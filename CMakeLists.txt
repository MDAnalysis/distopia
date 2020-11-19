cmake_minimum_required(VERSION 3.12)
project(DISTOPIA LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# compiler flags
set(CMAKE_CXX_FLAGS "-mavx -O3 -ffast-math")

# options
option(DISTOPIA_BUILD_TEST "Build distopia testing binary." ON)
option(DISTOPIA_BUILD_TIMINGS "Build distopia timings binary." ON)


# add the directory structure 
set(DISTOPIA_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
set(DISTOPIA_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src/lib/src)
set(DISTOPIA_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/src/lib/include)
set(DISTOPIA_TEST_DIR ${CMAKE_CURRENT_LIST_DIR}/src/lib/tests)
set(DISTOPIA_COMPARE_DIR ${CMAKE_CURRENT_LIST_DIR}/src/compare)
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})


file(GLOB distopia_src ${DISTOPIA_SRC_DIR}/*.cpp )
file(GLOB compare_src ${DISTOPIA_COMPARE_DIR}/*.cpp )

list (APPEND distopia_src ${compare_src})

# find the level of SSE instruction supported
find_package(SSE)
if(C_SSE1_FOUND)
  message(STATUS "SSE1 <xmmintrin.h> found ")
  set(CMAKE_C_FLAGS "${C_SSE1_FLAGS} -DUSE_SSE1 ${CMAKE_C_FLAGS}")
endif()
if(C_SSE2_FOUND)
  message(STATUS "SSE2 <emmintrin.h> found")
  set(CMAKE_C_FLAGS "${C_SSE2_FLAGS} -DUSE_SSE2 ${CMAKE_C_FLAGS}")
endif()
if(C_SSE3_FOUND)
  message(STATUS "SSE3 <pmmintrin.h> found")
  set(CMAKE_C_FLAGS "${C_SSE3_FLAGS} -DUSE_SSE3 ${CMAKE_C_FLAGS}")
endif()
if(C_SSSE3_FOUND)
  message(STATUS "SSSE3 <tmmintrin.h> found")
  set(CMAKE_C_FLAGS "${C_SSSE3_FLAGS} -DUSE_SSSE3 ${CMAKE_C_FLAGS}")
endif()
if(C_SSE4_1_FOUND)
  message(STATUS "SSE4_1 <smmintrin.h> found")
  set(CMAKE_C_FLAGS "${C_SSE4_1_FLAGS} -DUSE_SSE4_1 ${CMAKE_C_FLAGS}")
endif()
if(C_SSE4_2_FOUND)
  message(STATUS "SSE4_2 <nmmintrin.h> found")
  set(CMAKE_C_FLAGS "${C_SSE4_2_FLAGS} -DUSE_SSE4_2 ${CMAKE_C_FLAGS}")
endif()

find_package(AVX)


# ARM specific flags
find_package(ARM)
if(NEON_FOUND)
  message(STATUS "Neon found with compiler flag : -mfpu=neon -D__NEON__")
  set(CMAKE_C_FLAGS "-mfpu=neon -D__NEON__ ${CMAKE_C_FLAGS}")
endif()
if(CORTEXA8_FOUND)
  message(STATUS "Cortex-A8 Found with compiler flag : -mcpu=cortex-a8")
  set(CMAKE_C_FLAGS "-mcpu=cortex-a8 -fprefetch-loop-arrays ${CMAKE_C_FLAGS}")
endif()
if(CORTEXA9_FOUND)
  message(STATUS "Cortex-A9 Found with compiler flag : -mcpu=cortex-a9")
  set(CMAKE_C_FLAGS "-mcpu=cortex-a9 ${CMAKE_C_FLAGS}")
endif()


include_directories(distopia ${DISTOPIA_INCLUDE_DIR})
include_directories(distopia ${DISTOPIA_COMPARE_DIR})
add_library(distopia SHARED ${distopia_src})

if(DISTOPIA_BUILD_TEST)
    # add googletest
    Include(GoogleTest)
    add_subdirectory(${DISTOPIA_TEST_DIR}/googletest)
    enable_testing()
    include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
    # add test executable
    add_executable(tests ${DISTOPIA_TEST_DIR}/tests.cpp)
    target_link_libraries(tests distopia)
    target_link_libraries(tests gtest gtest_main)
    # discover tests so that make test will work
    gtest_discover_tests(tests WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

if(DISTOPIA_BUILD_TIMINGS)
    add_executable(timings ${DISTOPIA_TEST_DIR}/timings.cpp)
    target_link_libraries(timings distopia)
endif()