cmake_minimum_required(VERSION 3.13)

project(
    "RIDMME"
    VERSION 0.1.0
    LANGUAGES C CXX)

# Global CMake variables are set here
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_WARNINGS "Enable to add warnings to a target." ON)
option(ENABLE_WARNINGS_AS_ERRORS "Enable to treat warnings as errors." OFF)
option(ENABLE_CLANG_FORMAT "Enable to add clang-format." ON)


# Project/Library Names
set(EXECUTABLE_NAME "example")
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/Source/Public")

set(RIDMME_LIBRARY_NAME "RIDMME")

# CMAKE MODULES

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/)

include(Warnings)
include(Tools)

add_clang_format_target()

find_package(bias)

# SUB DIRECTORIES
add_subdirectory(Source)
#add_subdirectory(Example)

