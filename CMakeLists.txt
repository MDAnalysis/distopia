cmake_minimum_required(VERSION 3.24)
project(distopia2_the_highway_warrior)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)

set(HWY_ENABLE_TESTS OFF)
add_subdirectory("highway")


add_executable(demo "test.cpp")
target_sources(demo PRIVATE "library.h" "library.cpp")
target_link_libraries(demo PRIVATE hwy)
target_include_directories(demo PUBLIC "highway")
target_include_directories(demo PUBLIC  ${CMAKE_SOURCE_DIR})


add_library(calc_bonds_lib)
target_sources(calc_bonds_lib PRIVATE "calc_bonds.h" "calc_bonds.cpp")
target_link_libraries(calc_bonds_lib PRIVATE hwy)
target_include_directories(calc_bonds_lib PUBLIC "highway")
target_include_directories(calc_bonds_lib PUBLIC  ${CMAKE_SOURCE_DIR})

add_executable(run_calc_bonds)
target_sources(run_calc_bonds PRIVATE "run_calc_bonds.cpp")
target_link_libraries(run_calc_bonds PUBLIC calc_bonds_lib)
target_include_directories(run_calc_bonds PUBLIC  ${CMAKE_SOURCE_DIR})
