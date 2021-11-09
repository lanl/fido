add_executable(main
  main.cpp
  logging.cpp
  mpi_context.cpp
  distributed_runner.cpp
  windows.cpp)

target_link_libraries(main PRIVATE
  NLopt::nlopt
  fmt::fmt
  spdlog::spdlog
  shoccs::shoccs
  MPI::MPI_CXX)
