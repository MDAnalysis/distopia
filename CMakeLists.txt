cmake_minimum_required(VERSION 3.12)
project(DISTOPIA LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)




# options
option(DISTOPIA_BUILD_TEST "Build distopia testing binary." on)
option(DISTOPIA_BUILD_TIMINGS "Build distopia timings binary." on)


# add the directory structure 
set(DISTOPIA_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
set(DISTOPIA_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src/lib/src)
set(DISTOPIA_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/src/lib/include)
set(DISTOPIA_TEST_DIR ${CMAKE_CURRENT_LIST_DIR}/src/lib/tests)
set(DISTOPIA_COMPARE_DIR ${CMAKE_CURRENT_LIST_DIR}/src/compare)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})

# include some pretty colors
include(${CMAKE_CURRENT_LIST_DIR}/cmake/ColorOutput.cmake)


# user options for SIMD instruction levels
option(DISTOPIA_USE_SSE1    "Use SSE1 instructions" on)
option(DISTOPIA_USE_SSE2    "Use SSE2 instructions" on)
option(DISTOPIA_USE_SSE3    "Use SSE3 instructions" on)
option(DISTOPIA_USE_SSSE3   "Use SSSE3 instructions" on)
option(DISTOPIA_USE_SSE4_1  "Use SSE4.1 instructions" on)
option(DISTOPIA_USE_SSE4_2  "Use SSE4.2 instructions" on)



set(CMAKE_REQUIRED_QUIET true)
# find the level of SSE instruction supported
find_package(SSE)
# find the level of AVX supported
find_package(AVX)
set(CMAKE_REQUIRED_QUIET false)




#check the levels of SIMD supported
if(C_SSE1_FOUND AND CXX_SSE1_FOUND)
  message(STATUS "SSE1 SUPPORTED for C and CXX")
  if(DISTOPIA_USE_SSE1)
  message(STATUS "${Cyan}DISTOPIA SSE1 SELECTED${ColourReset}")
  endif()
endif()
if(C_SSE2_FOUND AND CXX_SSE2_FOUND)
  message(STATUS "SSE2 SUPPORTED for C and CXX")
  if(DISTOPIA_USE_SSE2)
  message(STATUS "${Cyan}DISTOPIA SSE2 SELECTED${ColourReset}")
  endif()
endif()
if(C_SSE3_FOUND AND CXX_SSE3_FOUND)
  message(STATUS "SSE3 SUPPORTED for C and CXX")
  if(DISTOPIA_USE_SSE3)
  message(STATUS "${Cyan}DISTOPIA SSE3 SELECTED${ColourReset}")
  endif()
endif()
if(C_SSSE3_FOUND AND CXX_SSSE3_FOUND)
  message(STATUS "SSSE3 SUPPORTED for C and CXX")
  if(DISTOPIA_USE_SSSE3)
  message(STATUS "${Cyan}DISTOPIA SSSE3 SELECTED${ColourReset}")
  endif()
endif()
if(C_SSE4_1_FOUND AND CXX_SSE4_1_FOUND)
  message(STATUS "SSE4_1 SUPPORTED for C and CXX")
  if(DISTOPIA_USE_SSE4_1)
  message(STATUS "${Cyan}DISTOPIA SSE4_1 SELECTED${ColourReset}")
  endif()
endif()
if(C_SSE4_2_FOUND AND CXX_SSE4_2_FOUND)
  message(STATUS "SSE4_2 SUPPORTED for C and CXX")
  if(DISTOPIA_USE_SSE4_2)
  message(STATUS "${Cyan}DISTOPIA SSE4_2 SELECTED${ColourReset}")
  endif()
endif()





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


#find all the files
file(GLOB distopia_src ${DISTOPIA_SRC_DIR}/*.cpp )
file(GLOB compare_src ${DISTOPIA_COMPARE_DIR}/*.cpp )
list (APPEND distopia_src ${compare_src})

# manually set compiler flags
set(CMAKE_CXX_FLAGS "-mavx -O3 -ffast-math")

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