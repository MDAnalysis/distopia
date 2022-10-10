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

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

# distopia c++ library
add_subdirectory(libdistopia)
# python layer
add_subdirectory(distopia)
