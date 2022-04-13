#include <cmath>
#include <legion.h>
#include <nlopt.hpp>

#include <vector>

#include <filesystem>
#include <fstream>
#include <sol/sol.hpp>

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

namespace fs = std::filesystem;

using namespace Legion;

Logger log_fido("fido");

enum TaskIDs {
    TOP_LEVEL_TASK_ID,
    TOP_LEVEL_NLOPT_TASK_ID,
};

void top_level_task(const Task* task,
                    const std::vector<PhysicalRegion>& regions,
                    Context ctx,
                    Runtime* runtime)
{
    log_fido.print("top_level...");
    // parse input and decide how many nlopt tasks to file up
    int n = 1;

    auto sz = fs::file_size("input.lua");
    auto total_sz = sz + sizeof(int) + 1;
    void* task_buf = new char[total_sz];
    char* file_buf = (char*)(task_buf) + sizeof(int);
    std::ifstream is{"input.lua"};
    is.read(file_buf, sz);
    file_buf[sz] = '\0';

    // fire off a bunch of top-level-nlopt tasks
    for (int i = 0; i < n; i++) {
        *((int*)task_buf) = i;
        TaskLauncher l(TOP_LEVEL_NLOPT_TASK_ID, TaskArgument(task_buf, total_sz));
        runtime->execute_task(ctx, l);
    }
}

struct my_constraint_data {
    double a, b;
};

double myfunc(unsigned n, const double* x, double* grad, void* my_func_data)
{
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / std::sqrt(x[1]);
    }

    return std::sqrt(x[1]);
}

double my_constraint(unsigned n, const double* x, double* grad, void* data)
{
    my_constraint_data& d = *(my_constraint_data*)data;
    auto&& [a, b] = d;

    if (grad) {
        grad[0] = 3 * a * (a * x[0] + b) * (a * x[0] + b);
        grad[1] = -1;
    }
    return (a * x[0] + b) * (a * x[0] + b) * (a * x[0] + b) - x[1];
}

void top_level_nlopt_task(const Task* task,
                          const std::vector<PhysicalRegion>& regions,
                          Context ctx,
                          Runtime* runtime)
{
    // assert(task->arglen == sizeof(int));
    int n = *(const int*)task->args;
    const char* buf = (const char*)task->args + sizeof(int);
    log_fido.print("top-level-nlopt task %d", n);

    sol::state lua{};
    lua.open_libraries(
        sol::lib::base, sol::lib::string, sol::lib::package, sol::lib::math);
    lua.safe_script(buf);

    // would this be the right place to "serialize" the input (maybe into a bunch of json
    // strings) which could be passed as an extra argument to the nlopt optimize function.
    // That function would then create the index spaces and do a bunch of index launches
    // for the simulation runs.  The results of the simulations would be processed via a
    // reduction task.  And those would be further processed via reduction
    // Each index space with a series of input parameters would also need a corresponding
    // output region.  There would need to be some kind of non-overlapping way of doing
    // the output region

    // initialize nlopt options

    std::vector<double> lb = {-HUGE_VAL, 0}; /* lower bounds */
    nlopt::opt opt{nlopt::LD_MMA, 2};
    opt.set_lower_bounds(lb);
    opt.set_min_objective(myfunc, NULL);
    std::vector<my_constraint_data> d{{2, 0}, {-1, 1}};
    opt.add_inequality_constraint(my_constraint, &d[0], 1.0e-8);
    opt.add_inequality_constraint(my_constraint, &d[1], 1.0e-8);
    opt.set_xtol_rel(1e-4);

    std::vector x = {1.234, 5.678};
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

    return Runtime::start(argc, argv);
}
