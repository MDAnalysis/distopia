cmake_minimum_required(VERSION 3.12)
project(Xdist)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-mavx -O3 -ffast-math")

add_executable(Xdist main.cpp)
