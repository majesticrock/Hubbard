# CMake for HubbardMeanField
# requires: MPI, Handler


cmake_minimum_required(VERSION 3.19)

project(Hubbard VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")

set(MKL_LINK static)
set(MKL_THREADING gnu_thread)
set(MKL_INTERFACE lp64)

find_package(MKL CONFIG QUIET HINTS $ENV{MKLROOT})
if (MKL_FOUND)
    set(USE_MKL ON CACHE BOOL "Use Intel MKL" FORCE)
else()
    set(USE_MKL OFF CACHE BOOL "Use Intel MKL" FORCE)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/MKLFlags.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/SelectClusterBuildOptions.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/../../PhdUtility/cmake/mrock-information.cmake)

include_directories(
    $ENV{HOME}/usr/local/include
)

# System libraries
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

find_package(mrock REQUIRED CONFIG)

# subdirs
add_subdirectory(sources/Hubbard/Helper)
add_subdirectory(sources/Handler)
add_subdirectory(sources/Hubbard/Models)
add_subdirectory(sources/Hubbard/DensityOfStates)

# Configure main
add_executable(HubbardMeanField 
    sources/HubbardMeanField.cpp
)

SET_COMPILER_FLAGS(HubbardMeanField)
target_link_libraries(HubbardMeanField PUBLIC
    Handler
    MPI::MPI_CXX
)

#target_link_options(HubbardMeanField PRIVATE "-Wl,--trace-symbol")