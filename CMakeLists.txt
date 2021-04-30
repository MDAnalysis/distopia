cmake_minimum_required(VERSION 3.12)
project(DISTOPIA LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_CXX_EXTENSIONS OFF)



# options
option(DISTOPIA_BUILD_TEST "Build distopia testing binary." on)
option(DISTOPIA_BUILD_TIMINGS "Build distopia timings binary." on)
option(DISTOPIA_COVERAGE "generate a code coverage report" off)
if(DISTOPIA_COVERAGE)
  set(CMAKE_BUILD_TYPE COVERAGE)
endif()

# add the directory structure 
set(DISTOPIA_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
set(DISTOPIA_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src/lib/src)
set(DISTOPIA_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/src/lib/include)
set(DISTOPIA_TEST_DIR ${CMAKE_CURRENT_LIST_DIR}/src/lib/tests)
set(DISTOPIA_COMPARE_DIR ${CMAKE_CURRENT_LIST_DIR}/src/compare)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})

# include some pretty colors
include(${CMAKE_CURRENT_LIST_DIR}/cmake/ColorOutput.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/CodeCoverage.cmake)


# user options for SIMD instructions
option(DISTOPIA_SELECT_SIMD "Let DISTOPIA choose instruction level" off)
option(DISTOPIA_USE_SSE1    "Use SSE1 instructions" off)
option(DISTOPIA_USE_SSE2    "Use SSE2 instructions" off)
option(DISTOPIA_USE_SSE3    "Use SSE3 instructions" off)
option(DISTOPIA_USE_SSSE3   "Use SSSE3 instructions" off)
option(DISTOPIA_USE_SSE4_1  "Use SSE4.1 instructions" off)
option(DISTOPIA_USE_SSE4_2  "Use SSE4.2 instructions" off)
option(DISTOPIA_USE_AVX     "Use AVX  instructions" off)
option(DISTOPIA_USE_AVX2    "Use AVX2  instructions" off)
option(DISTOPIA_USE_AVX512  "Use AVX512  instructions" off)



set(SSE_OPT_LIST "")
list(APPEND SSE_OPT_LIST ${DISTOPIA_USE_SSE1};${DISTOPIA_USE_SSE2};
${DISTOPIA_USE_SSE3};${DISTOPIA_USE_SSSE3};${DISTOPIA_USE_SSE4_1};${DISTOPIA_USE_SSE4_2};
${DISTOPIA_USE_AVX};${DISTOPIA_USE_AVX2};${DISTOPIA_USE_AVX512})

set(_nopt 0)
set(_find_simd_flavour false)
foreach(_option ${SSE_OPT_LIST})
  if(_option)
    math(EXPR _nopt "${_nopt}+1")
  endif()
endforeach()

if(DISTOPIA_SELECT_SIMD)
  if(${_nopt} EQUAL 1)
  message(FATAL_ERROR "${Red}More than one SIMD flavour specified.${ColourReset}")
  endif()
  message(STATUS "${Green}DISTOPIA_SELECT_SIMD is set to ON, DISTOPIA will attempt to choose SIMD flavour for you.${ColourReset}")
endif()

if(NOT DISTOPIA_SELECT_SIMD)
  if(${_nopt} EQUAL 0)
    message(STATUS "${Yellow}No instruction set selected and DISTOPIA_SELECT_SIMD is set to OFF, DISTOPIA will NOT attempt to choose for you.${ColourReset}")
  elseif(${_nopt} EQUAL 1)
    message(STATUS "${Green}Instruction set selected by user.${ColourReset}")
    message(STATUS "${Green}Setting DISTOPIA_SELECT_SIMD to OFF.${ColourReset}")
  endif()
endif()

message(STATUS "Checking what level of SIMD is supported by compiler")
#check the levels of SIMD supported
set(CMAKE_REQUIRED_QUIET true)
find_package(SSE)
find_package(AVX)
set(CMAKE_REQUIRED_QUIET false)


# check both CXX and C supported and append to list for heuristics
#TODO turn this disgusting copy paste job into a macro
set(HIGHEST_SIMD_SUPPORTED "None")

if(C_SSE1_FOUND AND CXX_SSE1_FOUND)
  message(STATUS "SSE1 SUPPORTED for C and CXX")
  set(SSE1_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE1")
  list(APPEND SSE_SUPPORTED_LIST  ${SSE1_SUPPORTED})
  add_compile_definitions(DISTOPIA_SSE1_SUPPORTED=1)
else()
  set(SSE1_SUPPORTED false)
endif()

if(C_SSE2_FOUND AND CXX_SSE2_FOUND)
  message(STATUS "SSE2 SUPPORTED for C and CXX")
  set(SSE2_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE2")
  list(APPEND SSE_SUPPORTED_LIST  ${SSE2_SUPPORTED})
  add_compile_definitions(DISTOPIA_SSE2_SUPPORTED=1)
else()
  set(SSE2_SUPPORTED false)
endif()

if(C_SSE3_FOUND AND CXX_SSE3_FOUND)
  message(STATUS "SSE3 SUPPORTED for C and CXX")
  set(SSE3_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE3")
  list(APPEND SSE_SUPPORTED_LIST  ${SSE3_SUPPORTED})
  add_compile_definitions(DISTOPIA_SSE3_SUPPORTED=1)
else()
  set(SSE3_SUPPORTED false)
endif()

if(C_SSSE3_FOUND AND CXX_SSSE3_FOUND)
  message(STATUS "SSSE3 SUPPORTED for C and CXX")
  set(SSSE3_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSSE3")
  list(APPEND SSE_SUPPORTED_LIST  ${SSSE3_SUPPORTED})
  add_compile_definitions(DISTOPIA_SSSE3_SUPPORTED=1)
else()
  set(SSSE3_SUPPORTED false)
endif()

if(C_SSE4_1_FOUND AND CXX_SSE4_1_FOUND)
  message(STATUS "SSE4_1 SUPPORTED for C and CXX")
  set(SSE4_1_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE4_1")
  list(APPEND SSE_SUPPORTED_LIST  ${SSE4_1_SUPPORTED})
  add_compile_definitions(DISTOPIA_SSE4_1_SUPPORTED=1)
else()
  set(SSE4_1_SUPPORTED false)
endif()

if(C_SSE4_2_FOUND AND CXX_SSE4_2_FOUND)
  message(STATUS "SSE4_2 SUPPORTED for C and CXX")
  set(SSE4_2_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE4_2")  
  list(APPEND SSE_SUPPORTED_LIST  ${SSE4_2_SUPPORTED})
  add_compile_definitions(DISTOPIA_SSE4_2_SUPPORTED=1)
else()
  set(SSE4_2_SUPPORTED false)
endif()

if(C_AVX_FOUND AND CXX_AVX_FOUND)
  message(STATUS "AVX SUPPORTED for C and CXX")
  set(AVX_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "AVX")  
  list(APPEND SSE_SUPPORTED_LIST  ${AVX_SUPPORTED})
  add_compile_definitions(DISTOPIA_AVX_SUPPORTED=1)
else()
  set(AVX_SUPPORTED false)
endif()

if(C_AVX2_FOUND AND CXX_AVX2_FOUND)
  message(STATUS "AVX2 SUPPORTED for C and CXX")
  set(AVX2_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "AVX2")  
  list(APPEND SSE_SUPPORTED_LIST  ${AVX2_SUPPORTED})
  add_compile_definitions(DISTOPIA_AVX2_SUPPORTED=1)
else()
  set(AVX_SUPPORTED false)
endif()


if(C_AVX512_FOUND AND CXX_AVX512_FOUND)
  message(STATUS "AVX512 SUPPORTED for C and CXX")
  set(AVX512_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "AVX512")  
  list(APPEND SSE_SUPPORTED_LIST  ${AVX512_SUPPORTED})
  add_compile_definitions(DISTOPIA_AVX512_SUPPORTED=1)
else()
  set(AVX512_SUPPORTED false)
endif()

message(STATUS "Highest level of SIMD supported is ${HIGHEST_SIMD_SUPPORTED}")

if(DISTOPIA_SELECT_SIMD)
  message(STATUS "Selecting SIMD")
  set(DISTOPIA_USE_${HIGHEST_SIMD_SUPPORTED} ON)  
endif()



# choose SIMD and set flags 
set(DISTOPIA_SIMD_FLAGS)
set(DISTOPIA_FALLTHROUGH false)
if(DISTOPIA_USE_SSE1 AND SSE1_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA SSE1 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SIMD_FLAGS "-msse")
  add_compile_definitions(DISTOPIA_USE_SSE1=1)

elseif(DISTOPIA_USE_SSE2 AND SSE2_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA SSE2 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SIMD_FLAGS "-msse2")
  add_compile_definitions(DISTOPIA_USE_SSE2=1)

elseif(DISTOPIA_USE_SSE3 AND SSE3_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA SSE3 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SIMD_FLAGS "-msse3")
  add_compile_definitions(DISTOPIA_USE_SSE3=1)

elseif(DISTOPIA_USE_SSSE3 AND SSSE3_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA SSSE3 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SIMD_FLAGS "-mssse3")
  add_compile_definitions(DISTOPIA_USE_SSE3=1)

elseif(DISTOPIA_USE_SSE4_1 AND SSE4_1_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA SSE4_1 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SIMD_FLAGS "-msse4.1")
  add_compile_definitions(DISTOPIA_USE_SSE4_1=1)

elseif(DISTOPIA_USE_SSE4_2 AND SSE4_2_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA SSE4_2 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SIMD_FLAGS "-msse4.2")
  add_compile_definitions(DISTOPIA_USE_SSE4_2=1)

elseif(DISTOPIA_USE_AVX AND AVX_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA AVX SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SIMD_FLAGS "-mavx")
  add_compile_definitions(DISTOPIA_USE_AVX=1)

elseif(DISTOPIA_USE_AVX2 AND AVX2_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA AVX2 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SIMD_FLAGS "-mavx2")
  add_compile_definitions(DISTOPIA_USE_AVX2=1)

elseif(DISTOPIA_USE_AVX512 AND AVX512_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA AVX512 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SIMD_FLAGS "-mavx512") # note that this is a bit platform specific
  add_compile_definitions(DISTOPIA_USE_AVX512=1)

else()
  set(DISTOPIA_FALLTHROUGH true)
  message(STATUS "${Cyan} No SIMD type explicitly selected, fallthrough to -O3 flags.${ColourReset}")
  add_compile_definitions(DISTOPIA_FALLTHROUGH=1)
endif()


#add more and get rid of them as we please
set(COMPILER_FLAGS)
set(COMPILER_FLAGS_DEBUG)
set(COMPILER_FLAGS_RELEASE)
set(COMPILER_FLAGS_COVERAGE)

# do cross platform compatibility here
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  if(DISTOPIA_FALLTHROUGH)
    list(APPEND CXX_FLAGS "-march=native" "-mtune=native")
  endif()
  list(APPEND CXX_FLAGS "-O3" "-mfma" "-fno-math-errno" "-fno-trapping-math" "-Wno-ignored-attributes")
  list(APPEND CXX_FLAGS_DEBUG  "-g" "-O0" "-Wall"  "-Wextra" )
  list(APPEND CXX_FLAGS_RELEASE "-Wno-unused")
  list(APPEND CXX_FLAGS_COVERAGE "-g" "-O0" "-Wall"  "-Wextra" "--coverage")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  if(DISTOPIA_FALLTHROUGH)
    list(APPEND CXX_FLAGS "-march=native" "-mtune=native")
  endif()
  list(APPEND CXX_FLAGS "-O3" "-mfma" "-fno-math-errno" "-fno-trapping-math")
  list(APPEND CXX_FLAGS_DEBUG   "-g" "-O0" "-Wall"  "-Wextra")
  list(APPEND CXX_FLAGS_RELEASE "-Wno-unused")
  list(APPEND CXX_FLAGS_COVERAGE "-g" "-O0" "-Wall"  "-Wextra" "--coverage")
else()
  message(FATAL_ERROR "${Red}Compiler not supported${ColourReset}")
endif()

set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS_DEBUG} --coverage")
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} --coverage")

list(APPEND CXX_FLAGS ${DISTOPIA_SIMD_FLAGS})
list(APPEND CXX_FLAGS_DEBUG ${DISTOPIA_SIMD_FLAGS})
list(APPEND CXX_FLAGS_RELEASE ${DISTOPIA_SIMD_FLAGS})
list(APPEND CXX_FLAGS_COVERAGE ${DISTOPIA_SIMD_FLAGS})


#find all the files
file(GLOB distopia_src ${DISTOPIA_SRC_DIR}/*.cpp )
file(GLOB compare_src ${DISTOPIA_COMPARE_DIR}/*.cpp )
list (APPEND distopia_src ${compare_src})


include_directories(distopia ${DISTOPIA_INCLUDE_DIR})
include_directories(distopia ${DISTOPIA_COMPARE_DIR})
add_library(distopia SHARED ${distopia_src})

target_compile_options(distopia
  PRIVATE
    ${CXX_FLAGS}
    "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:Coverage>:${CXX_FLAGS_COVERAGE}>"
  )

if(DISTOPIA_BUILD_TEST)
  # add googletest
  Include(GoogleTest)
  add_subdirectory(${DISTOPIA_TEST_DIR}/googletest)
  enable_testing()
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
  # add test executable for internals
  add_executable(tests ${DISTOPIA_TEST_DIR}/tests.cpp)
  target_compile_options(tests
    PRIVATE
      ${CXX_FLAGS}
      "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
      "$<$<CONFIG:Coverage>:${CXX_FLAGS_COVERAGE}>"
    )
  target_link_libraries(tests distopia)
  target_link_libraries(tests gtest gtest_main)
  # discover tests so that make test will work
  gtest_discover_tests(tests WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  # test executable for actual distances
  add_executable(test_kernels  ${DISTOPIA_TEST_DIR}/test_kernels.cpp)
  target_compile_options(test_kernels
    PRIVATE
      ${CXX_FLAGS}
      "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
      "$<$<CONFIG:Coverage>:${CXX_FLAGS_COVERAGE}>"
    )
  target_link_libraries(test_kernels distopia)
  target_link_libraries(test_kernels gtest gtest_main)
  gtest_discover_tests(test_kernels WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


  add_executable(test_mda_match  ${DISTOPIA_TEST_DIR}/test_mda_match.cpp)
  target_compile_options(test_mda_match
    PRIVATE
      ${CXX_FLAGS}
      "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
      "$<$<CONFIG:Coverage>:${CXX_FLAGS_COVERAGE}>"
    )
  target_link_libraries(test_mda_match distopia)
  target_link_libraries(test_mda_match gtest gtest_main)
  # discover tests so that make test will work
  gtest_discover_tests(test_mda_match WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

if(DISTOPIA_BUILD_TIMINGS)
  if(NOT DISTOPIA_BUILD_TEST)
    message(FATAL_ERROR "Building googlebench timings requires building googletest also")
  endif()

  # old timings executable
  add_subdirectory(${DISTOPIA_TEST_DIR}/googlebench)
  add_executable(timings ${DISTOPIA_TEST_DIR}/timings.cpp)
  target_compile_options(timings
    PRIVATE
      ${CXX_FLAGS}
      "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
      "$<$<CONFIG:Coverage>:${CXX_FLAGS_COVERAGE}>"
    )
  target_link_libraries(timings distopia)

  # googlebench executable
  add_executable(benchmarks ${DISTOPIA_TEST_DIR}/benchmark.cpp)
  target_compile_options(benchmarks
  PRIVATE
    ${CXX_FLAGS}
    "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:Coverage>:${CXX_FLAGS_COVERAGE}>"
  )
  target_link_libraries(benchmarks distopia)
  target_link_libraries(benchmarks benchmark::benchmark)
endif()

if(DISTOPIA_COVERAGE) 
  message(STATUS "excluding ${PROJECT_SOURCE_DIR}/src/lib/tests/googletest/* from coverage")
  message(STATUS "excluding ${PROJECT_SOURCE_DIR}/src/lib/tests/googlebench/* from coverage")

  setup_target_for_coverage_gcovr_html(
    NAME tests_coverage                    
    EXECUTABLE tests  
    DEPENDENCIES tests
    EXCLUDE "${PROJECT_SOURCE_DIR}/src/lib/tests/googletest/*" "${PROJECT_SOURCE_DIR}/src/lib/tests/googlebench/*")


  setup_target_for_coverage_gcovr_html(
    NAME test_kernels_coverage                    
    EXECUTABLE test_kernels  
    DEPENDENCIES test_kernels
    EXCLUDE "${PROJECT_SOURCE_DIR}/src/lib/tests/googletest/*" "${PROJECT_SOURCE_DIR}/src/lib/tests/googlebench/*")
endif()