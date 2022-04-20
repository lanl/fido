#include <cmath>

#include "driver.hpp"

#include <fmt/core.h>
#include <fmt/ranges.h>

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
    int n = 1;

    // something like:
    auto d = driver("input.lua");

    // fire off a bunch of top-level-nlopt tasks
    for (int i = 0; i < n; i++) {
        //*((int*)task_buf) = i;
        // can we instead do
        // TaskLauncher l(TOP_LEVEL_NLOPT_TASK_ID, interface);
        // where we have provided a conversion operator?
        TaskLauncher l(TOP_LEVEL_NLOPT_TASK_ID, d);
        runtime->execute_task(ctx, l);

        // or should this be an index launch?
    }
}

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

    log_fido.print(
        fmt::format("running objective with params: {}", fmt::join(params, ", "))
            .c_str());

    log_fido.print("task arglen %lu", task->arglen);

    // need to update dr with the alphas
    dr.set_data(params);

    // get the number of things to run and launch an index task
    Rect<1> launch_bounds(0, dr.simulation_size() - 1);

    ArgumentMap arg_map;
    for (int i = 0; i < dr.simulation_size(); i++) arg_map.set_point(i, dr);

    IndexTaskLauncher index_launcher(
        SIMULATION_TASK_ID, launch_bounds, TaskArgument(NULL, 0), arg_map);

    FutureMap fm = runtime->execute_index_space(ctx, index_launcher);
    fm.wait_all_results();

    std::vector<double> res(dr.simulation_size());
    for (int i = 0; i < dr.simulation_size(); i++)
        res[i] = fm.get_result<double>(i);

    return dr.result(res);
}

void top_level_nlopt_task(const Task* task,
                          const std::vector<PhysicalRegion>& regions,
                          Context ctx,
                          Runtime* runtime)
{
    // assert(task->arglen == sizeof(int));
    int n = *(const int*)task->args;
    log_fido.print("top-level-nlopt task %d", n);

    // here we would maybe do something like
    auto dr = driver::from_task(task);
    auto opt = dr.opt(log_fido);
    auto x = dr.guess();

    objective_data obj_d{dr, task, regions, ctx, runtime};

    opt.set_max_objective(objective, &obj_d);

    double maxval;
    auto opt_res = opt.optimize(x, maxval);

    log_fido.print("task %d found maximum in %d evaluations at f(%g, %g) = %g",
                   n,
                   opt.get_numevals(),
                   x[0],
                   x[1],
                   maxval);

    // call nlopt.optimize
    // thee objective we pass to nlopt will launch a bunch of other tasks
    // when we are done optimizing, recursively launch this task
}

double simulation_task(const Task* task,
                       const std::vector<PhysicalRegion>& regions,
                       Context ctx,
                       Runtime* runtime)
{
    auto dr = driver::from_task(task);
    auto x = dr.params();

    log_fido.print(fmt::format("running simulation task {}: with params {}",
                               task->index_point.point_data[0],
                               fmt::join(x, ", "))
                       .c_str());
    return dr.run(task->index_point.point_data[0]);
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
