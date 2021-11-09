#include "windows.hpp"
#include <mpi.h>

namespace fido
{

//
// simulation_info_window
//
simulation_info_window::simulation_info_window(const mpi_context& ctx, int n)
    : n{n}, root{ctx.root()}
{
    ctx.make_root_win(w, &base, n);
}

void simulation_info_window::set_info(std::span<const double> x)
{
    if (root) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, w);
        MPI_Put(&x[0], x.size(), MPI_DOUBLE, 0, 0, x.size(), MPI_DOUBLE, w);
        // for (int i = 0; i < n; i++) base[i] = x[i];
        //  std::copy_n(x.begin(), n, base);
        MPI_Win_unlock(0, w);
    }
}

void simulation_info_window::get_info(std::span<double> x)
{
    MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, w);
    MPI_Get(x.data(), n, MPI_DOUBLE, 0, 0, n, MPI_DOUBLE, w);
    MPI_Win_unlock(0, w);
}

//
// simulation_index_window
//
simulation_index_window::simulation_index_window(const mpi_context& ctx, int flag)
    : root{ctx.root()}, flag{flag}
{
    ctx.make_root_win(w, &base, 1, flag);
}

void simulation_index_window::set_initial_index(int idx)
{
    if (root) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, w);
        MPI_Put(&idx, 1, MPI_INT, 0, 0, 1, MPI_INT, w);
        MPI_Win_unlock(0, w);
    }
}

void simulation_index_window::reset()
{
    if (root) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, w);
        MPI_Put(&flag, 1, MPI_INT, 0, 0, 1, MPI_INT, w);
        MPI_Win_unlock(0, w);
    }
}

//
// result_acc_window
//
result_acc_window::result_acc_window(const mpi_context& ctx) : root{ctx.root()}
{
    ctx.make_root_win(w, &base, 1);
}

void result_acc_window::simulation_finished()
{
    int inc = 1;
    int ret = 0;
    MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, w);
    MPI_Fetch_and_op(&inc, &ret, MPI_INT, 0, 0, MPI_SUM, w);
    MPI_Win_unlock(0, w);
}

void result_acc_window::reset()
{
    if (root) {
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, w);
        int i = 0;
        MPI_Put(&i, 1, MPI_INT, 0, 0, 1, MPI_INT, w);
        MPI_Win_unlock(0, w);
    }
}

//
// result_window
//
result_window::result_window(const mpi_context& ctx, int n)
    : nresults{n}, root{ctx.root()}
{
    ctx.make_root_win(w, &base, nresults);
}

void result_window::set_result(double x, int index)
{
    MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, w);
    MPI_Put(&x, 1, MPI_DOUBLE, 0, index, 1, MPI_DOUBLE, w);
    MPI_Win_unlock(0, w);
}

std::span<const double> result_window::data() const
{
    return std::span<const double>(base, nresults);
}

} // namespace fido
