cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(test)
include_directories(../include)
add_library(quadsieve ../src/sieve.cpp ../src/num_theory)
add_executable(test test.cpp)
target_link_libraries(test quadsieve)
