# Hubbard/CMakeLists.txt
cmake_minimum_required(VERSION 3.30)

project(Hubbard VERSION 1.0)

#
# Set local settings and set the default build type as Release
#

list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/usr/local")
set(Boost_DIR /usr/lib64/openmpi/lib/cmake/Boost-1.78.0)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'.")

#
# Library specifics
#

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
