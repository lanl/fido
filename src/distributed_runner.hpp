#pragma once

#include <optional>
#include <numeric>

#include <sol/sol.hpp>

#include "windows.hpp"

    namespace fido
{
struct sim_result {
    int num_evals;
    double max_val;
    std::vector<double> x;
};

class distributed_runner
{
    const mpi_context* ctx;
    sol::table opt;
    sol::table sim;
    sol::table cons;
    bool root;
    simulation_info_window sim_info;
    simulation_index_window sim_index;
    result_acc_window n_results;
    result_window results;

public:
    distributed_runner(const mpi_context&,
                       sol::table opt,
                       sol::table sim,
                       sol::table cons);

    std::optional<sim_result> run();
};
} // namespace fido
