#include <cmath>

#include "driver.hpp"

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fstream>
#include <string_view>

/*
 * Some thoughts on moving our custom mpi/coroutine systems to a Legion task based
 * model The top level task would have to: parse input, spawn the tasks that will
 * initialize nlopt
 *
 * the nlopt tasks should spawn other tasks from calls to the objective function
 * nlopt would only need to be run as a single task but should then spawn other tasks
 * using the data from nlopt (i.e. alpha, )
 *
 * combining the results of the objective function evaluations should be a reduction
 * task
 *
 * Some thoughts input abstractions.
 *
 * 1. A main input serialize routine that can be used to launch multiple independent
 * NLOPT_TASKs
 * 2. Within a NLOPT_TASK, need to construct an nlopt::opt from the serialized input and
 * call opt.optimize
 * 3. That call will dispatch to an objective function which will need access to the
 * serialized input inorder to create create an index space corresponding to the different
 * classes of simulations.  An index lauch here with serialized input chunks corresponding
 * to the simulation classes.  From there, we create another index space corresponding to
 * the simulations to be carried out for each class of simulations. So we need something
 * like `serialize_whole_input`, `serialize_simulation_classes`, `serialize_simulation
 * input`
 *
 */

using namespace Legion;

Logger log_fido("fido");

enum TaskIDs {
    TOP_LEVEL_TASK_ID,
    TOP_LEVEL_NLOPT_TASK_ID,
    SIMULATION_TASK_ID,
};

void top_level_task(const Task* task,
                    const std::vector<PhysicalRegion>& regions,
                    Context ctx,
                    Runtime* runtime)
{
    log_fido.print("top_level...");

    // parse input and decide how many nlopt tasks to file up
    std::string input_file = "input.lua";
    int n = 1;

    // This is how the legion examples parse input... I hate this.
    const InputArgs& args = Runtime::get_input_args();
    for (int i = 1; i < args.argc; i++) {
        if (!strcmp(args.argv[i], "-n")) n = atoi(args.argv[++i]);
        if (args.argv[i][0] != '-') {
            std::string_view file_arg{args.argv[i]};
            if (file_arg.ends_with(std::string_view(".lua"))) input_file = file_arg;
        }
    }

    auto d = driver(input_file);

    // need to pass a unique task index for solution logging
    Rect<1> launch_bounds(0, n - 1);
    ArgumentMap arg_map;
    for (int j = 0; j < n; j++) {
        d.task_index() = j;
        arg_map.set_point(j, d);
    }
    IndexTaskLauncher index_launcher(
        TOP_LEVEL_NLOPT_TASK_ID, launch_bounds, TaskArgument(NULL, 0), arg_map);
    runtime->execute_index_space(ctx, index_launcher);
}

//
// Need to package up our input data and legion specific data so we can launch
// a bunch of simulation tasks from the objective function.  Note that the objective
// function will be called by indirectly via nlopt rather than directly by us
//
struct objective_data {
    driver_span& dr;
    const Task* task;
    const std::vector<PhysicalRegion>& regions;
    Context& ctx;
    Runtime* runtime;
};

double objective(unsigned n, const double* x, double* grad, void* data)
{
    auto& obj_d = *reinterpret_cast<objective_data*>(data);
    auto&& [dr, task, regions, ctx, runtime] = obj_d;

    std::span<const double> params{x, n};

    // need to update dr with the alphas
    dr.set_data(params);

    // get the outer number of things to run and launch an index task
    auto sz = dr.simulation_size();

    // launch all tasks and collect futures
    std::vector<FutureMap> fm;
    fm.reserve(sz);

    for (int i = 0; i < sz; i++) {
        Rect<1> launch_bounds(0, dr.simulation_size(i) - 1);

        // probably doesn't make much sense to do an argmap with the same data to every
        // task
        ArgumentMap arg_map;
        for (int j = 0; j < dr.simulation_size(i); j++) arg_map.set_point(j, dr);

        IndexTaskLauncher index_launcher(
            SIMULATION_TASK_ID, launch_bounds, TaskArgument(&i, sizeof(int)), arg_map);

        fm.push_back(runtime->execute_index_space(ctx, index_launcher));
    }

    // process all tasks
    std::vector<double> res;
    res.reserve(sz);
    for (int i = 0; i < sz; i++) {
        fm[i].wait_all_results();
        std::vector<double> local_res(dr.simulation_size(i));

        for (int j = 0; j < dr.simulation_size(i); j++)
            local_res[j] = fm[i].get_result<double>(j);

        res.push_back(dr.result(i, local_res));
    }

    double result = dr.result(res);
    log_fido.debug(fmt::format("objective with params: {}\n>>> result: {}",
                               fmt::join(params, ", "),
                               result)
                       .c_str());
    return result;
}

//
// Serial constraint function called directly from the top-level-nlopt task
//
double constraint(unsigned n, const double* x, double*, void* data)
{
    auto& obj_d = *reinterpret_cast<objective_data*>(data);
    auto& dr = obj_d.dr;

    std::span<const double> params{x, n};

    // need to update dr with the alphas
    dr.set_data(params);

    // do in serial for now since the eigenvalue calculations are fast
    double result = dr.constraint();
    return result;
}

//
//
//
void top_level_nlopt_task(const Task* task,
                          const std::vector<PhysicalRegion>& regions,
                          Context ctx,
                          Runtime* runtime)
{

    auto dr = driver::from_task(task);

    int n = dr.task_index();
    log_fido.debug("top-level-nlopt task %d", n);

    auto opt = dr.opt(log_fido);
    auto x = dr.guess();

    objective_data obj_d{dr, task, regions, ctx, runtime};

    opt.set_max_objective(objective, &obj_d);
    opt.add_inequality_constraint(constraint, &obj_d, 0.0);

    double maxval;
    // nlopt will call the objective and constraint functions
    opt.optimize(x, maxval);

    // log results to appropriate file
    std::string file = dr.accept(maxval) ? fmt::format("success.{:06d}", n)
                                         : fmt::format("fail.{:06d}", n);

    std::ofstream out{file, std::ios::binary | std::ios::app};
    out.write(reinterpret_cast<const char*>(x.data()), x.size() * sizeof(double));
    out.write(reinterpret_cast<const char*>(&maxval), sizeof(maxval));
    out.close();

    log_fido.debug(fmt::format("task {} found maximum in {} evaluations at\n f({}) = {}",
                               n,
                               opt.get_numevals(),
                               fmt::join(x, ", "),
                               maxval)
                       .c_str());

    // launch another task if time permits
    if (Realm::Clock::current_time() < dr.time_limit()) {
        dr.task_index() = n;
        TaskLauncher launcher(TOP_LEVEL_NLOPT_TASK_ID, dr);
        runtime->execute_task(ctx, launcher);
    }
}

//
// Leaf task which runs the simulation and returns the results to the nlopt driver
//
double simulation_task(const Task* task,
                       const std::vector<PhysicalRegion>& regions,
                       Context ctx,
                       Runtime* runtime)
{
    int idx;
    auto dr = driver::from_task(idx, task);
    return dr.run(idx, task->index_point.point_data[0]);
}

int main(int argc, char* argv[])
{
    Runtime::set_top_level_task_id(TOP_LEVEL_TASK_ID);
    {
        TaskVariantRegistrar r(TOP_LEVEL_TASK_ID, "top_level");
        r.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
        Runtime::preregister_task_variant<top_level_task>(r, "top_level");
    }
    {
        TaskVariantRegistrar r(TOP_LEVEL_NLOPT_TASK_ID, "top_level_nlopt");
        r.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
        Runtime::preregister_task_variant<top_level_nlopt_task>(r, "top_level_nlopt");
    }

    {
        TaskVariantRegistrar r(SIMULATION_TASK_ID, "simulation");
        r.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
        r.set_leaf();
        Runtime::preregister_task_variant<double, simulation_task>(r, "simulation");
    }

    return Runtime::start(argc, argv);
}
