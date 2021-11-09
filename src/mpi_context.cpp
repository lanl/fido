#include "mpi_context.hpp"

namespace fido
{

    mpi_context::mpi_context() : comm_free{false}
    {
        comm = MPI_COMM_WORLD;
        MPI_Comm_rank(comm, &rank_);
        root_ = rank_ == 0;
    }

    mpi_context::~mpi_context()
    {
        if (comm_free) MPI_Comm_free(&comm);
    }
}
