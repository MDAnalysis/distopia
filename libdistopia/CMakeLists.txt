# Distopia cmake file 

# options for which binaries to build 
option(DISTOPIA_BUILD_TEST "Build distopia testing binaries." off)
option(DISTOPIA_BUILD_TIMINGS "Build distopia timings binary." off)

# include some pretty colors
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})
include(${CMAKE_CURRENT_LIST_DIR}/cmake/ColorOutput.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/CodeCoverage.cmake)


# user options for SIMD instructions
option(DISTOPIA_AUTO_SELECT_SIMD "Let Distopia choose highest instruction set only" off)
option(DISTOPIA_MANUAL_SELECT_SIMD "Build Distopia for a single instruction set which must be explicitly specified  " off)
option(DISTOPIA_AGRESSIVE_MARCH "Allow distopia to build with march=native and mtune=native" off)

option(DISTOPIA_DISPATCH "Build Distopia with multiple instruction sets" off)
option(DISTOPIA_DISPATCH_MAX "Build Distopia with as many instruction sets as is supported" off)
option(DISTOPIA_DISPATCH_MANUAL "Explicitly specify which instruction sets to use in dispatch" off)


# options for explicit selection of SIMD specify which ones
option(DISTOPIA_USE_SSE1    "Use SSE1 instructions" off) # VCL2 instrset 1
option(DISTOPIA_USE_SSE2    "Use SSE2 instructions" off) # VCL2 instrset 2
option(DISTOPIA_USE_SSE3    "Use SSE3 instructions" off) # VCL2 instrset 3
option(DISTOPIA_USE_SSSE3   "Use SSSE3 instructions" off) # VCL2 instrset 4
option(DISTOPIA_USE_SSE4_1  "Use SSE4.1 instructions" off) # VCL2 instrset 5
option(DISTOPIA_USE_SSE4_2  "Use SSE4.2 instructions" off) # VCL2 instrset 6
option(DISTOPIA_USE_AVX     "Use AVX  instructions" off) # VCL2 instrset 7
option(DISTOPIA_USE_AVX2    "Use AVX2  instructions" off) # VCL2 instrset 8
option(DISTOPIA_USE_AVX512  "Use AVX512  instructions" off) # VCL2 instrset 9+

# find if any of the manual options were specified
set(SSE_OPT_LIST "")
list(APPEND SSE_OPT_LIST ${DISTOPIA_USE_SSE1};${DISTOPIA_USE_SSE2};
${DISTOPIA_USE_SSE3};${DISTOPIA_USE_SSSE3};${DISTOPIA_USE_SSE4_1};${DISTOPIA_USE_SSE4_2};
${DISTOPIA_USE_AVX};${DISTOPIA_USE_AVX2};${DISTOPIA_USE_AVX512})

set(_nopt 0)
foreach(_option ${SSE_OPT_LIST})
  if(_option)
    math(EXPR _nopt "${_nopt}+1")
  endif()
endforeach()

# default case
if(NOT DISTOPIA_MANUAL_SELECT_SIMD AND NOT DISTOPIA_DISPATCH AND NOT DISTOPIA_DISPATCH_MAX AND NOT DISTOPIA_DISPATCH_MANUAL)
  message(STATUS "${Yellow}No CMake build options specified, defaulting to DISTOPIA_AUTO_SELECT_SIMD ${ColourReset}")
  set(DISTOPIA_AUTO_SELECT_SIMD on)
endif()


# deal with options for DISTOPIA_AUTO_SELECT_SIMD
if(DISTOPIA_AUTO_SELECT_SIMD)
  message(STATUS "${Green}DISTOPIA_AUTO_SELECT_SIMD is set to ON, Distopia will attempt to choose highest SIMD flavour ONLY.${ColourReset}")
  if(DISTOPIA_DISPATCH OR DISTOPIA_DISPATCH_MAX OR DISTOPIA_DISPATCH_MANUAL)
    message(FATAL_ERROR "${Red}DISTOPIA_DISPATCH is incompatable with auto selecting highest SIMD flavour.${ColourReset}")
  endif()
  if(${_nopt} GREATER 0)
    message(FATAL_ERROR "${Red}DISTOPIA_AUTO_SELECT_SIMD but an instruction set was specified manually.${ColourReset}")
  endif()
  if(DISTOPIA_MANUAL_SELECT_SIMD)
    message(FATAL_ERROR "${Red}DISTOPIA_AUTO_SELECT_SIMD is incompatable with DISTOPIA_MANUAL_SELECT_SIMD ${ColourReset}")
  endif()
endif()

# deal with options for DISTOPIA_MANUAL_SELECT_SIMD
if(DISTOPIA_MANUAL_SELECT_SIMD)
  message(STATUS "${Green}DISTOPIA_MANUAL_SELECT_SIMD is set to ON, Distopia expects a single instruction set provided on command line.${ColourReset}")
  
  if(DISTOPIA_DISPATCH OR DISTOPIA_DISPATCH_MAX OR DISTOPIA_DISPATCH_MANUAL)
    message(FATAL_ERROR "${Red}DISTOPIA_DISPATCH is incompatable with manually selecting SIMD flavour.${ColourReset}")
  endif()  
  if(${_nopt} EQUAL 0)
  message(FATAL_ERROR "${Red}No SIMD flavour specified explicitly.${ColourReset}")
  endif()
  if(${_nopt} EQUAL 1)
    message(STATUS "${Green}A single instruction set has been provided ${ColourReset}")
  else()
    message(FATAL_ERROR "${Red}More than one SIMD flavour specified explicitly.${ColourReset}")
  endif()
endif()

if(NOT DISTOPIA_DISPATCH AND DISTOPIA_DISPATCH_MAX)
  message(STATUS "${Yellow}DISTOPIA_DISPATCH_MAX selected but DISTOPIA_DISPATCH not set, setting DISTOPIA_DISPATCH to ON${ColourReset}")
  set(DISTOPIA_DISPATCH on)
endif()

if(NOT DISTOPIA_DISPATCH AND DISTOPIA_DISPATCH_MANUAL)
  message(STATUS "${Yellow}DISTOPIA_DISPATCH_MANUAL selected but DISTOPIA_DISPATCH not set, setting DISTOPIA_DISPATCH to ON${ColourReset}")
  set(DISTOPIA_DISPATCH on)
endif()

# deal with options for DISTOPIA_DISPATCH
if(DISTOPIA_DISPATCH)
  message(STATUS "${Green}DISTOPIA_DISPATCH is set to ON, DISTOPIA will be compiled for multiple instruction sets ${ColourReset}")
  if(DISTOPIA_AUTO_SELECT_SIMD OR DISTOPIA_MANUAL_SELECT_SIMD)
    message(FATAL_ERROR "${Red}Building dispatch library is incompatable with automatically or manually selecting SIMD flavour.${ColourReset}")
  endif()

  if(DISTOPIA_DISPATCH_MAX AND DISTOPIA_DISPATCH_MANUAL)
  message(FATAL_ERROR "${Red}DISTOPIA_DISPATCH_MAX is incompatable with DISTOPIA_DISPATCH_MANUAL ${ColourReset}")
  endif()

  if(DISTOPIA_DISPATCH_MANUAL)
    message(STATUS "${Green}DISTOPIA_DISPATCH_MANUAL selected, Distopia expects dispatch instruction sets to be provided on the command line ${ColourReset}")
    if(${_nopt} EQUAL 0)
    message(FATAL_ERROR "${Red}No SIMD flavour specified explicitly.${ColourReset}")
    elseif(${_nopt} GREATER_EQUAL 1)
      message(STATUS "${Green}At least one one SIMD flavour specified explicitly.${ColourReset}")
    endif()
    message(STATUS "${Green}For manual dispatch, minimum instruction set (SSE1) must be active, setting DISTOPIA_USE_SSE1 to ON${ColourReset}")
    set(DISTOPIA_USE_SSE1 ON)
  elseif(DISTOPIA_DISPATCH_MAX)
    message(STATUS "${Green}DISTOPIA_DISPATCH_MAX selected, Distopia will be build with to dispatch with all available instruction sets ${ColourReset}")

  elseif(NOT DISTOPIA_DISPATCH_MANUAL AND NOT DISTOPIA_DISPATCH_MAX)
    message(STATUS "${Yellow}DISTOPIA_DISPATCH selected but MANUAL or MAX dispatch not set, setting DISTOPIA_DISPATCH_MAX to ON${ColourReset}")
    message(STATUS "${Green}DISTOPIA_DISPATCH_MAX selected, Distopia will be build with to dispatch with all available instruction sets ${ColourReset}")
    set(DISTOPIA_DISPATCH_MAX ON)
  endif()
endif()



message(STATUS "${Cyan}Checking what level of SIMD is supported on this machine${ColourReset}")
#check the levels of SIMD supported
find_package(SSE)
find_package(AVX)


# check both CXX and C supported and append to list for heuristics
#TODO turn this disgusting copy paste job into a macro
set(HIGHEST_SIMD_SUPPORTED "None")

if(C_SSE1_FOUND AND CXX_SSE1_FOUND)
  message(STATUS "SSE1 SUPPORTED for C and CXX")
  set(SSE1_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE1")
  list(APPEND SSE_SUPPORTED_LIST  ${SSE1_SUPPORTED})
else()
  set(SSE1_SUPPORTED false)
endif()

if(C_SSE2_FOUND AND CXX_SSE2_FOUND)
  message(STATUS "SSE2 SUPPORTED for C and CXX")
  set(SSE2_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE2")
  list(APPEND SSE_SUPPORTED_LIST  ${SSE2_SUPPORTED})
else()
  set(SSE2_SUPPORTED false)
endif()

if(C_SSE3_FOUND AND CXX_SSE3_FOUND)
  message(STATUS "SSE3 SUPPORTED for C and CXX")
  set(SSE3_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE3")
  list(APPEND SSE_SUPPORTED_LIST  ${SSE3_SUPPORTED})
else()
  set(SSE3_SUPPORTED false)
endif()

if(C_SSSE3_FOUND AND CXX_SSSE3_FOUND)
  message(STATUS "SSSE3 SUPPORTED for C and CXX")
  set(SSSE3_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSSE3")
  list(APPEND SSE_SUPPORTED_LIST  ${SSSE3_SUPPORTED})
else()
  set(SSSE3_SUPPORTED false)
endif()

if(C_SSE4_1_FOUND AND CXX_SSE4_1_FOUND)
  message(STATUS "SSE4_1 SUPPORTED for C and CXX")
  set(SSE4_1_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE4_1")
  list(APPEND SSE_SUPPORTED_LIST  ${SSE4_1_SUPPORTED})
else()
  set(SSE4_1_SUPPORTED false)
endif()

if(C_SSE4_2_FOUND AND CXX_SSE4_2_FOUND)
  message(STATUS "SSE4_2 SUPPORTED for C and CXX")
  set(SSE4_2_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "SSE4_2")  
  list(APPEND SSE_SUPPORTED_LIST  ${SSE4_2_SUPPORTED})
else()
  set(SSE4_2_SUPPORTED false)
endif()

if(C_AVX_FOUND AND CXX_AVX_FOUND)
  message(STATUS "AVX SUPPORTED for C and CXX")
  set(AVX_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "AVX")  
  list(APPEND SSE_SUPPORTED_LIST  ${AVX_SUPPORTED})
else()
  set(AVX_SUPPORTED false)
endif()

if(C_AVX2_FOUND AND CXX_AVX2_FOUND)
  message(STATUS "AVX2 SUPPORTED for C and CXX")
  set(AVX2_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "AVX2")  
  list(APPEND SSE_SUPPORTED_LIST  ${AVX2_SUPPORTED})
else()
  set(AVX2_SUPPORTED false)
endif()


if(C_AVX512_FOUND AND CXX_AVX512_FOUND)
  message(STATUS "AVX512 SUPPORTED for C and CXX")
  set(AVX512_SUPPORTED true)
  set(HIGHEST_SIMD_SUPPORTED "AVX512")  
  list(APPEND SSE_SUPPORTED_LIST  ${AVX512_SUPPORTED})
else()
  set(AVX512_SUPPORTED false)
endif()

message(STATUS "${Cyan}Highest level of SIMD supported is ${HIGHEST_SIMD_SUPPORTED} ${ColourReset}")


if(DISTOPIA_AUTO_SELECT_SIMD)
  message(STATUS "${Green}DISTOPIA_AUTO_SELECT_SIMD is ON, using ${HIGHEST_SIMD_SUPPORTED} ${ColourReset}")
  set(DISTOPIA_USE_${HIGHEST_SIMD_SUPPORTED} on)  
endif()


set(all_instr
SSE1  
SSE2  
SSE3  
SSSE3 
SSE4_1
SSE4_2
AVX   
AVX2  
AVX512
)

if(DISTOPIA_DISPATCH_MAX)
message(STATUS "${Green}DISTOPIA_DISPATCH_MAX selected, setting all instruction sets to ON ${ColourReset}")
foreach(opt ${all_instr})
  set(DISTOPIA_USE_${opt} ON)
endforeach()
endif()

# choose SIMD and set flags 
set(DISTOPIA_SIMD_FLAGS)

if(DISTOPIA_USE_SSE1)
  message(STATUS "${Cyan}DISTOPIA SSE1 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SSE1_SIMD_FLAGS  "-msse")
  list(APPEND DISTOPIA_SIMD_FLAGS ${DISTOPIA_SSE1_SIMD_FLAGS})
endif()

if(DISTOPIA_USE_SSE2)
  message(STATUS "${Cyan}DISTOPIA SSE2 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SSE2_SIMD_FLAGS "-msse2")
  list(APPEND DISTOPIA_SIMD_FLAGS ${DISTOPIA_SSE2_SIMD_FLAGS})
endif()

if(DISTOPIA_USE_SSE3)
  message(STATUS "${Cyan}DISTOPIA SSE3 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SSE3_SIMD_FLAGS "-msse3")
  list(APPEND DISTOPIA_SIMD_FLAGS ${DISTOPIA_SSE3_SIMD_FLAGS})
endif()

if(DISTOPIA_USE_SSSE3)
  message(STATUS "${Cyan}DISTOPIA SSSE3 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SSSE3_SIMD_FLAGS "-mssse3")
  list(APPEND DISTOPIA_SIMD_FLAGS ${DISTOPIA_SSSE3_SIMD_FLAGS})
endif()

if(DISTOPIA_USE_SSE4_1 AND SSE4_1_SUPPORTED)
  message(STATUS "${Cyan}DISTOPIA SSE4_1 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SSE4_1_SIMD_FLAGS "-msse4.1")
  list(APPEND DISTOPIA_SIMD_FLAGS ${DISTOPIA_SSE4_1_SIMD_FLAGS})
endif()

if(DISTOPIA_USE_SSE4_2)
  message(STATUS "${Cyan}DISTOPIA SSE4_2 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_SSE4_2_SIMD_FLAGS "-msse4.2")
  list(APPEND DISTOPIA_SIMD_FLAGS ${DISTOPIA_SSE4_2_SIMD_FLAGS})
endif()

if(DISTOPIA_USE_AVX)
  message(STATUS "${Cyan}DISTOPIA AVX SELECTED${ColourReset}")
  list(APPEND DISTOPIA_AVX_SIMD_FLAGS "-mavx")
  list(APPEND DISTOPIA_SIMD_FLAGS ${DISTOPIA_AVX_SIMD_FLAGS})
endif()

if(DISTOPIA_USE_AVX2)
  message(STATUS "${Cyan}DISTOPIA AVX2 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_AVX2_SIMD_FLAGS "-mavx2"  "-mfma")
  list(APPEND DISTOPIA_SIMD_FLAGS ${DISTOPIA_AVX2_SIMD_FLAGS})
endif()

if(DISTOPIA_USE_AVX512)
  message(STATUS "${Cyan}DISTOPIA AVX512 SELECTED${ColourReset}")
  list(APPEND DISTOPIA_AVX512_SIMD_FLAGS "-mavx512f" "-mfma") # note that this is a bit platform specific
  list(APPEND DISTOPIA_SIMD_FLAGS ${DISTOPIA_AVX512_SIMD_FLAGS}) # note that this is a bit platform specific
endif()

if(NOT DISTOPIA_DISPATCH)
message(STATUS "${Green}DISTOPIA_DISPATCH is OFF. Are we using march and mtune flags?${ColourReset}")
  if(DISTOPIA_AGRESSIVE_MARCH)
    message(STATUS "${Green}DISTOPIA_AGRESSIVE_MARCH is ON. Enabling aggressive optimisations for CPU type with -march=native and -mtune=native ${ColourReset}")
    message(STATUS "${Yellow}Warning this will OVERRIDE SIMD specific options ${ColourReset}")
    list(APPEND CXX_FLAGS "-march=native" "-mtune=native" "-mfma")
    message(STATUS "${Cyan}SIMD FLAGS are ${CXX_FLAGS} ${ColourReset}")
  else()
    message(STATUS "${Green}Aggressive CPU type optimisations (march=native mtune=native) disabled, enable them with DISTOPIA_AGRESSIVE_MARCH${ColourReset}")
    list(APPEND CXX_FLAGS ${DISTOPIA_SIMD_FLAGS})
    list(APPEND CXX_FLAGS_DEBUG ${DISTOPIA_SIMD_FLAGS})
    list(APPEND CXX_FLAGS_RELEASE ${DISTOPIA_SIMD_FLAGS})
    message(STATUS "${Cyan}SIMD FLAGS are ${CXX_FLAGS}${ColourReset}")
  endif()
endif()


# do cross platform compatibility here
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  list(APPEND CXX_FLAGS "-O3" "-fno-math-errno" "-fno-trapping-math" "-Wno-ignored-attributes")
  list(APPEND CXX_FLAGS_DEBUG  "-g" "-O0" "-Wall"  "-Wextra" )
  list(APPEND CXX_FLAGS_RELEASE "-Wno-unused")
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  list(APPEND CXX_FLAGS "-O3" "-fno-math-errno" "-fno-trapping-math")
  list(APPEND CXX_FLAGS_DEBUG   "-g" "-O0" "-Wall"  "-Wextra")
  list(APPEND CXX_FLAGS_RELEASE "-Wno-unused")
else()
  message(FATAL_ERROR "${Red}Compiler not supported${ColourReset}")
endif()


# add the directory structure 
set(DISTOPIA_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
set(DISTOPIA_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/libdistopia/lib)
set(DISTOPIA_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/libdistopia/lib)
set(DISTOPIA_HEADER_DIR ${CMAKE_CURRENT_LIST_DIR}/libdistopia/include)
set(VCL2_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/vcl2)
set(DISTOPIA_TEST_DIR ${CMAKE_CURRENT_LIST_DIR}/tests)
set(DISTOPIA_COMPARE_DIR ${CMAKE_CURRENT_LIST_DIR}/compare)


#find all the files
file(GLOB distopia_src ${DISTOPIA_SRC_DIR}/*.cpp )
file(GLOB vcl_src ${VCL2_INCLUDE_DIR}/instrset_detect.cpp)



include_directories(libdistopia ${DISTOPIA_INCLUDE_DIR})
include_directories(libdistopia ${DISTOPIA_COMPARE_DIR})
include_directories(libdistopia ${VCL2_INCLUDE_DIR})

# we just build the library once with all source and our set of SIMD FLAGS
if(NOT DISTOPIA_DISPATCH)
  add_library(libdistopia SHARED ${distopia_src})

# need an object library for each ISA
elseif(DISTOPIA_DISPATCH)
  message(STATUS "${Magenta}Creating targets for dispatch${ColourReset}")
  # add the main library
  add_library(libdistopia SHARED)

  foreach(opt ${all_instr})
    if(DISTOPIA_USE_${opt})
      message(STATUS "${Magenta}Creating dispatch target for ${opt} target name distopia_${opt} with flags ${DISTOPIA_${opt}_SIMD_FLAGS} ${ColourReset}")
      if(${opt} STREQUAL SSE1)
      # a bunch of special things need to be embedded in the mandatory SSE1 build. 
      # 1. we need to compile the sources for instrset_detect.cpp so that we can use it
      # 2. we need to make available a bunch of flags that tell us what the available SIMD flags are
       list(APPEND  distopia_src ${vcl_src})
        message(STATUS "${Green}Appending nessecary VCL sources for target = distopia_${opt}${ColourReset}")
      endif()
      add_library(distopia_${opt} OBJECT ${distopia_src})
      target_compile_options(distopia_${opt} PRIVATE ${DISTOPIA_${opt}_SIMD_FLAGS})
      target_compile_definitions(distopia_${opt} PRIVATE DISTOPIA_USE_${opt}=1)
      target_compile_definitions(distopia_${opt} PUBLIC DISTOPIA_DISPATCH=1)
      if(${opt} STREQUAL SSE1)
      list(REMOVE_ITEM  distopia_src ${vcl_src})
      message(STATUS "${Green}Removing nessecary VCL sources for target > distopia_${opt}${ColourReset}")
      foreach(opt_internal ${all_instr})
      if(DISTOPIA_USE_${opt_internal})
        target_compile_definitions(distopia_SSE1 PRIVATE DISTOPIA_${opt_internal}_AVAILABLE=1)
      endif()
      endforeach()
      endif()
      target_link_libraries(libdistopia distopia_${opt})
    else()
      message(STATUS "${Yellow}Dispatch inactive for ${opt} ${ColourReset}")
    endif()
  endforeach()



else()
# should never happen but just in case 
message(FATAL_ERROR "Invalid option combination specified")
endif()

target_compile_options(libdistopia
  PRIVATE
    ${CXX_FLAGS}
    "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
  )

if(DISTOPIA_BUILD_TEST)
  # add googletest
  Include(GoogleTest)
  add_subdirectory(${DISTOPIA_TEST_DIR}/googletest)
  enable_testing()
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
  # add test executable for internals
  file(GLOB tests_src ${DISTOPIA_COMPARE_DIR}/*.cpp )
  list (APPEND tests_src ${DISTOPIA_TEST_DIR}/unit_tests.cpp)
  add_executable(unit_tests ${tests_src})
  target_include_directories(unit_tests PRIVATE ${DISTOPIA_SOURCE_DIR})
  target_compile_options(unit_tests
    PRIVATE
      ${CXX_FLAGS}
      "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
    )
  target_link_libraries(unit_tests libdistopia)
  target_link_libraries(unit_tests gtest gtest_main)

  # discover tests so that make test will work
  gtest_discover_tests(unit_tests WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

  add_executable(test_matches  ${DISTOPIA_TEST_DIR}/test_matches.cpp)
  target_compile_options(test_matches
    PRIVATE
      ${CXX_FLAGS}
      "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
      "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
    )
  target_link_libraries(test_matches libdistopia)
  target_link_libraries(test_matches gtest gtest_main)
  # discover tests so that make test will work
  gtest_discover_tests(test_matches WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

endif()

if(DISTOPIA_BUILD_TIMINGS)
  if(NOT DISTOPIA_BUILD_TEST)
    message(FATAL_ERROR "Building googlebench timings requires building googletest also")
  endif()

  SET(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Turn off googlebench tests" FORCE)
  # google bench
  add_subdirectory(${DISTOPIA_TEST_DIR}/googlebench)

  # googlebench executable
  file(GLOB benchmarks_src ${DISTOPIA_COMPARE_DIR}/*.cpp )
  list (APPEND benchmarks_src ${DISTOPIA_TEST_DIR}/benchmark.cpp)
  add_executable(benchmarks ${benchmarks_src}) 
  list(APPEND CXX_FLAGS "-msse3")
  target_compile_options(benchmarks
  PRIVATE
    ${CXX_FLAGS}
    "$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>"
  )
  target_link_libraries(benchmarks libdistopia)
  target_link_libraries(benchmarks benchmark::benchmark)
endif()

set_target_properties(libdistopia
        PROPERTIES PUBLIC_HEADER ${DISTOPIA_ROOT_SOURCE_DIR}/libdistopia/include/distopia.h)

# cmake will append the name lib already , avoid liblibdistopia
set_target_properties(libdistopia PROPERTIES OUTPUT_NAME distopia)

install(TARGETS libdistopia
        DESTINATION distopia)
