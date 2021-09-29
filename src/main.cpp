#include <cmath>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <functional>
#include <limits>
#include <nlopt.hpp>
#include <shoccs.hpp>
#include <sol/sol.hpp>
#include <span>
#include <vector>

#include <cxxopts.hpp>

#include <spdlog/spdlog.h>

/*
  Run a list of constraints from a sol::table with the following format

  {
    -- array of all the simulations to run as part of this constraint
        simulations = {},
    -- sets appropriate values in simulation table index
        set_values = function (span<double>) end,
    -- returns the appropriate real value from the simulation result
        result = function (real3) return real; end
    -- aggregrates all the results
        aggregate = function (span<double>) return real; end
}
*/
struct constraint_runner {
    sol::table sims;
    std::vector<double> results;
    sol::function set_values;
    sol::function result;
    sol::function agg;

    constraint_runner() = default;
    constraint_runner(const sol::table& t)
        : sims{t["simulations"]},
          results(sims.size()),
          set_values{t["set_values"]},
          result{t["result"]},
          agg{t["aggregate"]}
    {
    }

    double operator()(std::span<const double> p)
    {
        for (int i = 0; i < (int)results.size(); i++) {
            set_values(i + 1, p);
            auto r = ccs::simulation_run(sims[i + 1]);
            assert(r);
            results[i] = result(*r);
        }

        return agg(results);
    }
};

double objective(unsigned n, const double* x, double*, void* data)
{
    sol::table& tbl = *reinterpret_cast<sol::table*>(data);

    {
        auto rng = std::span<const double>(x, n);
        spdlog::info("running objective with params: {}", fmt::join(rng, ", "));
    }

    double max_time = tbl["step_controller"]["max_time"];
    for (unsigned i = 0; i < n; i++) { tbl["scheme"]["alpha"][i + 1] = x[i]; }

    auto result = ccs::simulation_run(tbl);
    assert(result);

    auto&& [time, error, _] = *result;

    tbl["constraint"] = max_time - time;

    if (time < max_time) {
        return 10 * time / max_time;
    } else {
        return 10 - std::log(error);
    }
}

double constraint(unsigned n, const double* x, double*, void* data)
{
    constraint_runner& runner = *reinterpret_cast<constraint_runner*>(data);

    std::span<const double> params{x, n};

    return runner(params);
}

int main(int argc, char* argv[])
{

    cxxopts::Options options("fido", "Run the finite-difference optimizer a given input");

    // clang-format off
    options.add_options()
        ("input-file", "Main lua input file", cxxopts::value<std::string>())
        ("help", "Print usage");
    // clang-format on
    options.parse_positional("input-file");

    auto result = options.parse(argc, argv);

    if (result.count("help") || result.arguments().size() == 0) {
        std::cout << options.help() << '\n';
        return 0;
    }

    sol::state lua;
    lua.open_libraries(
        sol::lib::base, sol::lib::string, sol::lib::package, sol::lib::math);

    if (!result.count("input-file")) {
        spdlog::error("input file must be specified");
        return 1;
    }

    lua.script_file(result["input-file"].as<std::string>());

    if (!lua["Simulation"].valid()) {
        spdlog::error("top level `Simulation` table must be specified");
        return 1;
    }
    sol::table tbl = lua["Simulation"];

    nlopt::opt opt(nlopt::LN_COBYLA, 9);

    // std::vector<double> lower_bounds{-5, 0};
    // std::vector<double> upper_bounds{5, 10};

    // opt.set_lower_bounds(lower_bounds);
    // opt.set_upper_bounds(upper_bounds);

    // int count{};
    // auto my_l_func = [&count](unsigned n, const double* x, double* grad) {
    //     ++count;
    //     if (grad) {
    //         grad[0] = 0.0;
    //         grad[1] = 0.5 / std::sqrt(x[1]);
    //     }
    //     return std::sqrt(x[1]);
    // };

    opt.set_max_objective(objective, &tbl);

    if (!lua["Constraints"][1].valid()) {
        spdlog::error("top level `Constraints` table must be specified");
        return 1;
    }

    sol::table t = lua["Constraints"][1];
    constraint_runner runner{t};
    opt.add_inequality_constraint(constraint, &runner, 0.0);
    // opt.add_inequality_constraint(my_constraint, &cd[1], 1e-8);

    opt.set_xtol_rel(1e-5);
    opt.set_xtol_abs(1e-8);
    opt.set_maxeval(50);
    opt.set_initial_step(0.1);

    auto x = std::vector(9, 0.0);
    double maxval;

    auto opt_res = opt.optimize(x, maxval);
    fmt::print("found maximum in {} evaluations at f({}) = {}\n",
               opt.get_numevals(),
               fmt::join(x, ", "),
               maxval);
}
