#pragma once

#include "mpi_context.hpp"
#include <cppcoro/generator.hpp>
#include <mpi.h>
#include <span>

namespace fido
{

// Wrapper type around an MPI_Win for the simulation info (currenty the parameters being
// optimized).  The root process will set the data and then set the initial index in the
// index window.  All process will then consume the parameters and run the simulations
// with them
class simulation_info_window
{
    MPI_Win w;
    double* base;
    int n;
    bool root;

public:
    simulation_info_window(const mpi_context&, int n);

    void set_info(std::span<const double> x);

    void get_info(std::span<double> x);

    int size() const { return n; }

    ~simulation_info_window() { MPI_Win_free(&w); }
};

// Wrapper type around an MPI_Win for the simulation index.  The root process is expected
// to call `set_initial_index` once the simulation info has been properly set.  All
// processes are to consume the indices via the generator interface via
// `simulation_indices`
class simulation_index_window
{

    MPI_Win w;
    int* base;
    bool root;
    int flag;

public:
    simulation_index_window(const mpi_context&, int flag = -1);

    void set_initial_index(int idx = 0);
    void reset();

    cppcoro::generator<int> simulation_indices(int n) const
    {
        // first check that the flag value has been overwritten
        int idx = flag;
        while (idx == flag) {
            MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, w);
            MPI_Get(&idx, 1, MPI_INT, 0, 0, 1, MPI_INT, w);
            MPI_Win_unlock(0, w);
        }

        // Loop through all valid indices
        int inc = 1;
        while (idx < n) {
            MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, w);
            MPI_Fetch_and_op(&inc, &idx, MPI_INT, 0, 0, MPI_SUM, w);
            MPI_Win_unlock(0, w);

            if (idx < n) co_yield idx;
        }
    }

    ~simulation_index_window() { MPI_Win_free(&w); }
};

// Wrapper type around an MPI_Win that keeps track of the total number of results computed
// The root process running nlopt will use this to decide when to return control to nlopt
class result_acc_window
{
    MPI_Win w;
    int* base;
    bool root;

public:
    result_acc_window(const mpi_context&);

    void simulation_finished();
    void reset();

    ~result_acc_window() { MPI_Win_free(&w); }

    int count() const
    {
        int res;
        MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, w);
        MPI_Get(&res, 1, MPI_INT, 0, 0, 1, MPI_INT, w);
        MPI_Win_unlock(0, w);
        return res;
    }
};

// Wrapper around an MPI_Win that contains a list of all results.  Each process will write
// their results to the corresponding simulation index
class result_window
{
    MPI_Win w;
    double* base;
    int nresults;
    bool root;

public:
    result_window(const mpi_context&, int n);

    void set_result(double, int index);

    int size() const { return nresults; }

    std::span<const double> data() const;

    ~result_window() { MPI_Win_free(&w); }
};
} // namespace fido
