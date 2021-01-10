cmake_minimum_required(VERSION 3.9)

project("one-solver" CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)

set(CMAKE_CXX_FLAGS "-fsycl")

find_package(Boost 1.65)

add_subdirectory(src)
add_subdirectory(tests)