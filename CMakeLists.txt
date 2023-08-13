cmake_minimum_required(VERSION 3.24)
project(distopia2_the_highway_warrior)

set(HWY "/home/richard/code/highway/build/libhwy.a")

add_executable(demo "test.cpp")
target_sources(demo PRIVATE "library.h" "library.cpp")
target_link_libraries(demo PRIVATE ${HWY})
target_include_directories(demo PUBLIC "/home/richard/code/highway/")


add_library(calc_bonds_lib)
target_sources(calc_bonds_lib PRIVATE "calc_bonds.h" "calc_bonds.cpp")
target_link_libraries(calc_bonds_lib PRIVATE ${HWY})
target_include_directories(calc_bonds_lib PUBLIC "/home/richard/code/highway/")

add_executable(run_calc_bonds)
target_sources(run_calc_bonds PRIVATE "run_calc_bonds.cpp")
target_link_libraries(run_calc_bonds PRIVATE calc_bonds_lib)