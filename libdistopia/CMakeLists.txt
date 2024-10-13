
set(HWY_ENABLE_TESTS OFF CACHE BOOL "Enable HWY tests" FORCE)


set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable googlebench tests" FORCE)
add_subdirectory("googlebench")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)


add_library(libdistopia)
target_sources(libdistopia PRIVATE "src/distopia.cpp")
# needed for M_PI symbol on Win
target_compile_definitions(libdistopia PUBLIC -D_USE_MATH_DEFINES)
target_include_directories(libdistopia PUBLIC
        "highway"
        ${CMAKE_CURRENT_LIST_DIR}
        "${CMAKE_CURRENT_LIST_DIR}/include")


# the link target for HWY differs if it is an external library
if(DISTOPIA_USE_EXTERNAL_HWY)
    message(STATUS "Using external HWY")
    find_package(HWY REQUIRED)
    target_link_libraries(libdistopia PRIVATE hwy::hwy)
else()
    message(STATUS "Using vendored HWY subpackage")
    add_subdirectory("highway")
    target_link_libraries(libdistopia PRIVATE hwy)
endif()
    

add_executable(bench)
target_sources(bench PRIVATE "test/bench.cpp")
target_link_libraries(bench PUBLIC benchmark::benchmark)
target_link_libraries(bench PUBLIC libdistopia)
target_include_directories(bench PUBLIC  ${CMAKE_CURRENT_LIST_DIR})


Include(GoogleTest)
add_subdirectory("googletest")
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
include_directories(${gmock_SOURCE_DIR}/include ${gmock_SOURCE_DIR})


add_executable(test)
target_sources(test PRIVATE "test/test.cpp")
target_link_libraries(test PUBLIC gtest gtest_main)
target_link_libraries(test PUBLIC gmock gmock_main)

target_link_libraries(test PUBLIC libdistopia)
target_include_directories(test PUBLIC  ${CMAKE_SOURCE_DIR})
gtest_discover_tests(test WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


add_executable(targets)
target_sources(targets PRIVATE "test/targets.cpp")
target_link_libraries(targets PUBLIC libdistopia)
target_include_directories(targets PUBLIC  ${CMAKE_SOURCE_DIR})

add_executable(test_mda_match)
target_sources(test_mda_match PRIVATE "test/test_mda_match.cpp")
target_link_libraries(test_mda_match PUBLIC gtest gtest_main)
target_link_libraries(test PUBLIC gmock gmock_main)
target_link_libraries(test_mda_match PUBLIC libdistopia)
target_include_directories(test_mda_match PUBLIC  ${CMAKE_SOURCE_DIR})
gtest_discover_tests(test_mda_match WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


# cmake will append the name lib already , avoid liblibdistopia
set_target_properties(libdistopia PROPERTIES OUTPUT_NAME distopia)
if(APPLE)
    set_target_properties(libdistopia PROPERTIES INSTALL_RPATH "@loader_path")
else()
    set_target_properties(libdistopia PROPERTIES INSTALL_RPATH "\$ORIGIN")
endif()
install(TARGETS libdistopia DESTINATION distopia)