#pragma once

#include "logging.hpp"
#include <algorithm>
#include <mpi.h>
#include <type_traits>

namespace fido
{

template <typename T>
concept Bytes = std::is_trivially_copyable_v<std::decay_t<T>>;

class mpi_context
{
    bool comm_free;
    bool root_;
    int rank_;

    MPI_Comm comm;

public:
    mpi_context();

    constexpr bool root() const { return root_; }
    constexpr bool rank() const { return rank_; }

    logs root_logger(const std::string& name) const { return logs{root(), name}; }

    template <Bytes T>
    int Bcast(T& t) const
    {
        return MPI_Bcast(reinterpret_cast<void*>(&t), sizeof(T), MPI_CHAR, 0, comm);
    }

    template <typename T>
    int Bcast(T* t, int sz) const
    {
        return MPI_Bcast(reinterpret_cast<void*>(t), sz * sizeof(*t), MPI_CHAR, 0, comm);
    }

    int barrier() const { return MPI_Barrier(comm); }

    template <typename T>
    void make_root_win(MPI_Win& w, T** base, int elements, T init = {}) const
    {
        auto sz = sizeof(T);
        if (root()) {
            MPI_Win_allocate(sz * elements, sz, MPI_INFO_NULL, comm, base, &w);
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, w);
            std::fill_n(*base, elements, init);
            MPI_Win_unlock(0, w);
        } else {
            MPI_Win_allocate(0, 1, MPI_INFO_NULL, comm, base, &w);
        }
    }

    ~mpi_context();
};
} // namespace fido
