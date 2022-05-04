add_executable(fido-exe
  main.cpp)

target_link_libraries(fido-exe
  NLopt::nlopt
  fmt::fmt
  Legion::Legion
  sol2::sol2
  shoccs::shoccs
  lua)

set_target_properties(fido-exe PROPERTIES OUTPUT_NAME "fido")
install(TARGETS fido-exe)
