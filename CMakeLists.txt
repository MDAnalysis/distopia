cmake_minimum_required(VERSION 3.22)
project(DISTOPIA LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)
include_directories(${NumPy_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/libdistopia/libdistopia/include)
include_directories(${PROJECT_SOURCE_DIR}/libdistopia/libdistopia/lib)
include_directories(${PROJECT_SOURCE_DIR}/libdistopia/vcl2)


set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}")

# distopia c++ library
add_subdirectory(libdistopia)
# python layer
add_subdirectory(distopia)
