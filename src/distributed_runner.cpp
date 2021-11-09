#include "distributed_runner.hpp"

#include <shoccs.hpp>

#include <nlopt.hpp>

namespace fido
{
/*
  Run a list of constraints from a sol::table with the following format

  {
    -- array of all the simulations to run as part of this constraint
        simulations = {},
    -- sets appropriate values in simulation table index
        set_values = function (self, i, span<double>) end,
    -- returns the appropriate real value from the simulation result
        result = function (self, real3) return real; end
    -- aggregrates all the results
        aggregate = function (self, span<double>) return real; end
}
*/
struct constraint_runner {
    sol::table cons;
    std::vector<double> results;
    sol::function set_values;
    sol::function result;
    sol::function agg;

    constraint_runner() = default;
    constraint_runner(sol::table t)
        : cons{t},
          results(static_cast<sol::table>(t["simulations"]).size()),
          set_values{cons["set_values"]},
          result{cons["result"]},
          agg{cons["aggregate"]}
    {
    }

    double operator()(std::span<const double> p)
    {
        for (int i = 0; i < (int)results.size(); i++) {
            set_values(cons, i + 1, p);
            auto r = ccs::simulation_run(cons["simulations"][i + 1]);
            assert(r);
            results[i] = result(cons, *r);
        }

        return agg(cons, results);
    }
};

double constraint(unsigned n, const double* x, double*, void* data)
{
    constraint_runner& runner = *reinterpret_cast<constraint_runner*>(data);

    std::span<const double> params{x, n};

    return runner(params);
}

/*
  Run a list of constraints from a sol::table with the following format

  {
    -- array of all the simulations to run as part of this constraint
        simulations = {},
    -- sets appropriate values in simulation table index
        set_values = function (self, i, span<double>) end,
    -- returns the appropriate real value from the simulation result
        result = function (self, real3) return real; end
    -- aggregrates all the results
        aggregate = function (self, span<double>) return real; end
}
*/
struct objective_runner {
    const mpi_context* ctx;
    sol::table sims;
    sol::function set_values;
    sol::function result;
    sol::function agg;
    simulation_info_window& info;
    simulation_index_window& s_i;
    result_acc_window& result_count;
    result_window& results;

    objective_runner(const mpi_context* ctx,
                     sol::table t,
                     simulation_info_window& info,
                     simulation_index_window& s_i,
                     result_acc_window& result_count,
                     result_window& results)
        : ctx{ctx},
          sims{t},
          set_values{sims["set_values"]},
          result{sims["result"]},
          agg{sims["aggregate"]},
          info{info},
          s_i{s_i},
          result_count{result_count},
          results{results}
    {
    }

    void operator()()
    {
        std::vector<double> p(info.size());

        bool run;
        ctx->Bcast(run);
        while (run) {

            for (auto&& i : s_i.simulation_indices(results.size())) {

                info.get_info(p);
                set_values(sims, i + 1, p);
                auto r = ccs::simulation_run(sims["simulations"][i + 1]);
                assert(r);
                results.set_result(result(sims, *r), i);
                result_count.simulation_finished();
            }
            // sync for reset
            ctx->barrier();
            // check if root is done
            ctx->Bcast(run);
        }
    }

    double operator()(std::span<const double> p)
    {
        // initial setup for x and index;
        // spdlog::info("setting info...");
        info.set_info(p);
        // spdlog::info("setting initial index...");
        s_i.set_initial_index(0);
        // tell waiting procs to run
        bool run = true;
        // spdlog::info("broadcasting run");
        ctx->Bcast(run);

        for (auto&& i : s_i.simulation_indices(results.size())) {
            // spdlog::info("i: {}\n", i);
            set_values(sims, i + 1, p);
            auto r = ccs::simulation_run(sims["simulations"][i + 1]);
            assert(r);
            // spdlog::info("setting result...");
            results.set_result(result(sims, *r), i);
            // spdlog::info("marking finished...");
            result_count.simulation_finished();
        }

        // need to wait for count to finish - MPI_Barrier?
        // spdlog::info("waiting...");
        ctx->barrier();

        auto x = results.data();
        // spdlog::info("results: {}", fmt::join(x, ", "));
        s_i.reset();
        result_count.reset();

        return agg(sims, x);
    }
};

double objective(unsigned n, const double* x, double*, void* data)
{
    objective_runner& runner = *reinterpret_cast<objective_runner*>(data);
    std::span<const double> params{x, n};
    spdlog::info("running objective with params: {}", fmt::join(params, ", "));

    return runner(params);
}

distributed_runner::distributed_runner(const mpi_context& ctx,
                                       sol::table opt,
                                       sol::table sim,
                                       sol::table cons)
    : ctx{&ctx},
      opt{opt},
      sim{sim},
      cons{cons},
      root{ctx.root()},
      sim_info{ctx, opt["dims"].get<int>()},
      sim_index{ctx},
      n_results{ctx},
      results(ctx, static_cast<sol::table>(sim["simulations"]).size())
{
}

std::optional<sim_result> distributed_runner::run()
{
    constraint_runner con_r{cons};
    objective_runner obj_r{ctx, sim, sim_info, sim_index, n_results, results};

    if (root) {

        // prepare nlopt objective
        nlopt::opt opt(nlopt::LN_COBYLA, sim_info.size());

        opt.set_max_objective(objective, &obj_r);

        opt.add_inequality_constraint(constraint, &con_r, 0.0);

        opt.set_xtol_rel(1e-5);
        opt.set_xtol_abs(1e-8);
        opt.set_maxeval(10);
        opt.set_initial_step(0.1);

        auto x = std::vector(sim_info.size(), 0.0);
        double maxval;

        auto opt_res = opt.optimize(x, maxval);

        // let the other procs know we are done
        bool run = false;
        ctx->Bcast(run);
        // logger(spdlog::level::info,
        //        "found maximum in {} evaluations at f({}) = {}\n",
        //        opt.get_numevals(),
        //        fmt::join(x, ", "),
        //        maxval);
        // add constraints
        return sim_result{opt.get_numevals(), maxval, x};
    } else {
        obj_r();
        return sim_result{};
    }
}
} // namespace fido
