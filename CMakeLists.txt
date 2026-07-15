# CMake for HubbardMeanField
# requires: MPI, Handler
cmake_minimum_required(VERSION 3.30)

project(Hubbard VERSION 1.0)

list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/usr/local")

# System libraries
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

find_package(mrock 
    REQUIRED 
    COMPONENTS 
        symbolic_operators 
        iEoM 
        utility
)

include(${CMAKE_CURRENT_LIST_DIR}/cmake/hubbard_options.cmake)

# subdirs
add_subdirectory(sources/Hubbard/Helper)
add_subdirectory(sources/Handler)
add_subdirectory(sources/Hubbard/Models)
add_subdirectory(sources/Hubbard/DensityOfStates)

# Configure main
add_executable(HubbardMeanField 
    sources/HubbardMeanField.cpp
)

target_link_libraries(HubbardMeanField PUBLIC
    Handler
    MPI::MPI_CXX
    ${mrock_LIBRARIES}
    hubbard_options
)

# provided by the mrock library
mrock_generate_information_header(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/mrock/info.h"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)