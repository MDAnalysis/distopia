
set(HWY_ENABLE_TESTS OFF CACHE BOOL "Enable HWY tests" FORCE)
add_subdirectory("highway")

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable googlebench tests" FORCE)
add_subdirectory("googlebench")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)


add_library(distopia)
target_sources(distopia PRIVATE "src/distopia.cpp")
target_link_libraries(distopia PRIVATE hwy)
target_include_directories(distopia PUBLIC
        "highway"
        ${CMAKE_CURRENT_LIST_DIR}
        "${CMAKE_CURRENT_LIST_DIR}/include")


add_executable(run_bench)
target_sources(run_bench PRIVATE "bench.cpp")
target_link_libraries(run_bench PUBLIC benchmark::benchmark)
target_link_libraries(run_bench PUBLIC distopia)
target_include_directories(run_bench PUBLIC  ${CMAKE_CURRENT_LIST_DIR})


Include(GoogleTest)
add_subdirectory("googletest")
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(testing)
target_sources(testing PRIVATE "testing.cpp")
target_link_libraries(testing PUBLIC gtest gtest_main)
target_link_libraries(testing PUBLIC distopia)
target_include_directories(testing PUBLIC  ${CMAKE_SOURCE_DIR})
