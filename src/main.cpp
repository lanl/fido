#include <cmath>

#include <cstdio>
#include <filesystem>

#include <functional>
#include <limits>

#include <sol/sol.hpp>
#include <span>
#include <spdlog/common.h>
#include <string>
#include <vector>

#include <cxxopts.hpp>

#include <mpi.h>

#include "distributed_runner.hpp"
#include "mpi_context.hpp"
#include "mpi_singleton_wrapper.hpp"

namespace fs = std::filesystem;

/* Process the input file by having the root process slurp the file and broadcast the
 * contents to the other processes.  Returns the lua state resulting from running the
 * input file */
sol::state slurp_input_file(const fido::mpi_context& ctx, std::string filename)
{

    int sz;
    if (ctx.root()) sz = fs::file_size(filename);
    ctx.Bcast(sz);
    char* buf = new char[sz + 1];

    if (ctx.root()) {
        FILE* f = fopen(filename.c_str(), "r");
        fread(buf, 1, sz, f);
        fclose(f);
        buf[sz] = 0;
    }
    ctx.Bcast(buf, sz + 1);
    std::string context{std::move(buf)};

    sol::state lua;
    lua.open_libraries(
        sol::lib::base, sol::lib::string, sol::lib::package, sol::lib::math);
    lua.script(context);

    assert(lua["NLopt"].valid());
    assert(lua["Constraints"].valid());
    assert(lua["Simulations"].valid());

    return lua;
}

// double objective(unsigned n, const double* x, double*, void* data)
// {
//     sol::table& tbl = *reinterpret_cast<sol::table*>(data);

//     {
//         auto rng = std::span<const double>(x, n);
//         spdlog::info("running objective with params: {}", fmt::join(rng, ", "));
//     }

//     double max_time = tbl["step_controller"]["max_time"];
//     for (unsigned i = 0; i < n; i++) { tbl["scheme"]["alpha"][i + 1] = x[i]; }

//     auto result = ccs::simulation_run(tbl);
//     assert(result);

//     auto&& [time, error, _] = *result;

//     tbl["constraint"] = max_time - time;

//     if (time < max_time) {
//         return 10 * time / max_time;
//     } else {
//         return 10 - std::log(error);
//     }
// }

int main(int argc, char* argv[])
{
    fido::mpi_global_env env{argc, argv};

    cxxopts::Options options("fido", "Run the finite-difference optimizer a given input");

    // clang-format off
    options.add_options()
        ("input-file", "Main lua input file", cxxopts::value<std::string>())
        ("help", "Print usage");
    // clang-format on
    options.parse_positional("input-file");

    auto result = options.parse(argc, argv);

    fido::mpi_context ctx{};
    auto logger = ctx.root_logger("root");

    if (result.count("help") || result.arguments().size() == 0) {
        logger(spdlog::level::info, "\n{}\n", options.help());
        return 0;
    }

    if (!result.count("input-file")) {
        logger(spdlog::level::err, "input file must be specified");
        return 1;
    }

    sol::state lua = slurp_input_file(ctx, result["input-file"].as<std::string>());

    if (!lua["Simulations"].valid()) {
        logger(spdlog::level::err, "top level `Simulations` table must be specified");
        return 1;
    }

    if (sol::table t = lua["Simulations"]; t.size() == 0) {
        logger(spdlog::level::err, "top level `Simulations` must not be empty");
        return 1;
    }

    if (!lua["Constraints"].valid()) {
        logger(spdlog::level::err, "top level `Constraints` table must be specified");
        return 1;
    }

    if (sol::table t = lua["Constraints"]; t.size() == 0) {
        logger(spdlog::level::err, "top level `Constraints` table must not be empty");
        return 1;
    }

    // hardcode one instance of Simulations/Constraints for now.  Need to generalize to
    // more
    sol::table opt = lua["NLopt"];
    sol::table sim = lua["Simulations"][1];
    sol::table cons = lua["Constraints"][1];
    {
        int dims = opt["dims"];
        logger(spdlog::level::info, "dims = {}\n", dims);
    }
    fido::distributed_runner dr{ctx, opt, sim, cons};

    auto res = dr.run();
    if (!res) {
        logger(spdlog::level::err, "runner failed");
        return 1;
    }

    auto [num_evals, max_val, x] = *res;

    logger(spdlog::level::info,
           "found maximum in {} evaluations at f({}) = {}\n",
           num_evals,
           fmt::join(x, ", "),
           max_val);

    // if (ctx.root()) {
    //     sol::table tbl = lua["Simulation"];

    //     nlopt::opt opt(nlopt::LN_COBYLA, 9);

    //     opt.set_max_objective(objective, &tbl);

    //     sol::table t = lua["Constraints"][1];
    //     constraint_runner runner{t};
    //     opt.add_inequality_constraint(constraint, &runner, 0.0);

    //     opt.set_xtol_rel(1e-5);
    //     opt.set_xtol_abs(1e-8);
    //     opt.set_maxeval(50);
    //     opt.set_initial_step(0.1);

    //     auto x = std::vector(9, 0.0);
    //     double maxval;

    //     auto opt_res = opt.optimize(x, maxval);
    //     logger(spdlog::level::info,
    //            "found maximum in {} evaluations at f({}) = {}\n",
    //            opt.get_numevals(),
    //            fmt::join(x, ", "),
    //            maxval);
    // } else {
    // }
}
