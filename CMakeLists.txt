cmake_minimum_required(VERSION 3.9)

project("one-solver" CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(ENV{I_MPI_CXX} "dpcpp")
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpiicpc)

find_package(Boost 1.65)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsycl -std=c++17")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lsycl -lboost_serialization")

add_subdirectory(src)
add_subdirectory(tests)