cmake_minimum_required(VERSION 3.16)

project(fido LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(NLopt REQUIRED)
find_package(fmt REQUIRED)

include(GNUInstallDirs)


add_subdirectory(docs)
add_subdirectory(src)
