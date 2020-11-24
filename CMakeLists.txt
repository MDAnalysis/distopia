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


# user options for SIMD instructions
option(DISTOPIA_SELECT_SIMD "Let DISTOPIA choose instruction level" off)
option(DISTOPIA_USE_SSE1    "Use SSE1 instructions" off)
option(DISTOPIA_USE_SSE2    "Use SSE2 instructions" off)
option(DISTOPIA_USE_SSE3    "Use SSE3 instructions" off)
option(DISTOPIA_USE_SSSE3   "Use SSSE3 instructions" off)
option(DISTOPIA_USE_SSE4_1  "Use SSE4.1 instructions" off)
option(DISTOPIA_USE_SSE4_2  "Use SSE4.2 instructions" off)


set(SSE_LEVEL_LIST "")
list(APPEND SSE_LEVEL_LIST DISTOPIA_USE_SSE1;DISTOPIA_USE_SSE2
DISTOPIA_USE_SSE3;DISTOPIA_USE_SSSE3};DISTOPIA_USE_SSE4_1;DISTOPIA_USE_SSE4_2)

set(SSE_OPT_LIST "")
list(APPEND SSE_OPT_LIST ${DISTOPIA_SELECT_SIMD};${DISTOPIA_USE_SSE1};${DISTOPIA_USE_SSE2};
${DISTOPIA_USE_SSE3};${DISTOPIA_USE_SSSE3};${DISTOPIA_USE_SSE4_1};${DISTOPIA_USE_SSE4_2})

set(_nopt 0)
set(_find_simd_flavour false)
foreach(_option ${SSE_OPT_LIST})
  if(_option)
    math(EXPR _nopt "${_nopt}+1")
  endif()
endforeach()

if(${_nopt} EQUAL 0)
  message(STATUS "${Green}No instruction set selected DISTOPIA will attempt to choose for you.${ColourReset}")
  message(STATUS "${Green}Setting DISTOPIA_SELECT_SIMD to ON.${ColourReset}")
  set(DISTOPIA_SELECT_SIMD on)
elseif(${_nopt} EQUAL 1)
  message(STATUS "${Green}Instruction set selected by user.${ColourReset}")
  message(STATUS "${Green}Setting DISTOPIA_SELECT_SIMD to OFF.${ColourReset}")
  set(DISTOPIA_SELECT_SIMD off)
else()
  message(FATAL_ERROR "${Red}More than one SIMD flavour specified.${ColourReset}")
endif()

#check the levels of SIMD supported
set(CMAKE_REQUIRED_QUIET true)
find_package(SSE)
find_package(AVX)
set(CMAKE_REQUIRED_QUIET false)

# check both CXX and C supported and append to list for heuristics
set(SSE_SUPPORTED_LIST "")

if(C_SSE1_FOUND AND CXX_SSE1_FOUND)
  message(STATUS "SSE1 SUPPORTED for C and CXX")
  set(SSE1_SUPPORTED true)
  if(DISTOPIA_USE_SSE1)
    message(STATUS "${Cyan}DISTOPIA SSE1 SELECTED${ColourReset}")
  else()
    list(APPEND SSE_SUPPORTED_LIST  ${SSE1_SUPPORTED})
  endif()
endif()

if(C_SSE2_FOUND AND CXX_SSE2_FOUND)
  message(STATUS "SSE2 SUPPORTED for C and CXX")
  set(SSE1_SUPPORTED true)
  if(DISTOPIA_USE_SSE2)
    message(STATUS "${Cyan}DISTOPIA SSE2 SELECTED${ColourReset}")
  else()
    list(APPEND SSE_SUPPORTED_LIST  ${SSE1_SUPPORTED})
  endif()
endif()

if(DISTOPA_DISTOPIA_SELECT_SIMD)
  foreach(_sse_level ${SSE_SUPPORTED_LIST})
  

endforeach()




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