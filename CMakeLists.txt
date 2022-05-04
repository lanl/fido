cmake_minimum_required(VERSION 3.16)

project(fido LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()


set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Lua REQUIRED)
# make lua a target
add_library(lua INTERFACE IMPORTED)
target_include_directories(lua INTERFACE ${LUA_INCLUDE_DIR})
target_link_libraries(lua INTERFACE ${LUA_LIBRARIES})

find_package(NLopt REQUIRED)
find_package(fmt 8 REQUIRED)
find_package(Legion REQUIRED)

find_package(shoccs REQUIRED)
find_package(sol2 REQUIRED)

include(GNUInstallDirs)

if(NOT APPLE)
    set(CMAKE_INSTALL_RPATH $ORIGIN)
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_subdirectory(src)
