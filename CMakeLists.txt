cmake_minimum_required(VERSION 3.16)

project(fido LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Lua REQUIRED)
# make lua a target
add_library(lua INTERFACE IMPORTED)
target_include_directories(lua INTERFACE ${LUA_INCLUDE_DIR})
target_link_libraries(lua INTERFACE ${LUA_LIBRARIES})

find_package(NLopt REQUIRED)
find_package(fmt 8 REQUIRED)
find_package(shoccs REQUIRED)
find_package(cppcoro REQUIRED)

# include(CMakeFindDependencyMacro)
# find_dependency(shoccs)
find_package(sol2 REQUIRED)
find_package(range-v3 REQUIRED)
find_package(pugixml 1.10 REQUIRED)
find_package(spdlog 1.9 REQUIRED)
find_package(cxxopts REQUIRED)
find_package(Boost REQUIRED) # for header only mp11
find_package(MKL REQUIRED)
find_package(MPI REQUIRED)

include(GNUInstallDirs)


add_subdirectory(docs)
add_subdirectory(src)
