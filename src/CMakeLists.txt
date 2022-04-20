add_executable(main
  main.cpp)

target_link_libraries(main PRIVATE
  NLopt::nlopt
  fmt::fmt
  Legion::Legion
  sol2::sol2
  shoccs::shoccs
  lua)
