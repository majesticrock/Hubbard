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

# provided by the mrock library
mrock_generate_information_header(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/hubbard/info.h"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUT_HEADER MROCK_GENERATED_INFO_HEADER
    OUT_INCLUDE_DIR MROCK_GENERATED_INCLUDE_DIR
)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/hubbard_options.cmake)
# Sets include dir for the generated metadata header file
target_include_directories(hubbard_options INTERFACE ${MROCK_GENERATED_INCLUDE_DIR})

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
