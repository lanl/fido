#pragma once

#include <mpi.h>

namespace fido {

// global guard for MPI environment ensuring MPI_Finalize gets called
struct mpi_global_env {
    mpi_global_env(int argc, char* argv[]) {
        MPI_Init(&argc, &argv);
    }
    ~mpi_global_env(){
        MPI_Finalize();
    }
};
}
