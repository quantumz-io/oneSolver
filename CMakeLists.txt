cmake_minimum_required(VERSION 3.9)

project("one-solver" CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(ENV{I_MPI_CXX} "dpcpp")
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpiicpc)

# For beta10 and later due to change in location for dpl library.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsycl -std=c++17 -I $ENV{ONEAPI_ROOT}/dpl/2021.1-beta10/linux/include -fsycl-unnamed-lambda")

# Beta 09
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fsycl -std=c++17 -fsycl-unnamed-lambda")


set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lsycl -lboost_serialization")


find_package(Boost 1.65)

add_subdirectory(src)
add_subdirectory(tests)