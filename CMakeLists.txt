# CMake for HubbardMeanField
# requires: MPI, Handler


cmake_minimum_required(VERSION 3.19)

project(HubbardMeanField VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")

if(CLUSTER_BUILD)
    include(${CMAKE_SOURCE_DIR}/cmake/ClusterCompilerFlags.cmake)
else()
    include(${CMAKE_SOURCE_DIR}/cmake/DefaultCompilerFlags.cmake)
endif()

include_directories(
    $ENV{HOME}/usr/local/include
)

# System libraries
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

# SymbolicOperators library
find_package(SymbolicOperators REQUIRED CONFIG)
# Utility library
find_package(Utility REQUIRED)


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