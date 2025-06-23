# CMake for HubbardMeanField
# requires: MPI, Handler


cmake_minimum_required(VERSION 3.19)

project(Hubbard VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")

if(CLUSTER_BUILD)
    #include(${CMAKE_SOURCE_DIR}/cmake/CascadelakeCompilerFlags.cmake)
    include(${CMAKE_SOURCE_DIR}/cmake/IcelakeCompilerFlags.cmake)
else()
    include(${CMAKE_SOURCE_DIR}/cmake/DefaultCompilerFlags.cmake)
endif()
include(${CMAKE_CURRENT_SOURCE_DIR}/../../PhdUtility/cmake/mrock-information.cmake)

include_directories(
    $ENV{HOME}/usr/local/include
)

option(USE_MKL "Use Intel MKL for LAPACK" OFF)
if(USE_MKL)
    list(APPEND CMAKE_PREFIX_PATH "/opt/intel/oneapi/mkl/latest")
    find_package(MKL CONFIG REQUIRED)
endif()

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