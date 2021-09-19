cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

project(tests)
include_directories(../include)
enable_testing()

add_executable(sieve_tests sieve_tests.cpp ../src/sieve.cpp ../src/num_theory.cpp)
target_link_libraries(
                    sieve_tests
                    gtest_main
                  )
include(GoogleTest)
gtest_discover_tests(sieve_tests)
