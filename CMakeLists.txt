cmake_minimum_required(VERSION 3.12)
project(distopia)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-mavx -O3 -ffast-math")

add_library(distopia SHARED distopia.cpp distopia.h vanilla.cpp vanilla.h calc_distances.h distancekernels.h)
add_executable(timings main.cpp)
target_link_libraries(timings distopia)