# HubbardMeanField/CMakeLists.txt

cmake_minimum_required(VERSION 3.19)

project(HubbardMeanField VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")

include_directories(
    $ENV{HOME}/usr/local/include
)

add_executable(HubbardMeanField 
    sources/HubbardMeanField.cpp
)

if(NOT WIN32)
    set(BoldCyan "${Esc}[1;36m")
    set(ColourReset "${Esc}[m")
endif()
if(CLUSTER_BUILD)
    include(${CMAKE_SOURCE_DIR}/cmake/ClusterCompilerFlags.cmake)
    message("${BoldCyan}Building for cascadelake!${ColourReset}")
else()
    include(${CMAKE_SOURCE_DIR}/cmake/DefaultCompilerFlags.cmake)
    message("${BoldCyan}Building for the local machine!${ColourReset}")
endif()
SET_COMPILER_FLAGS(HubbardMeanField)

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS serialization iostreams)
find_package(ZLIB REQUIRED)
find_package(Utility REQUIRED)
find_package(SymbolicOperators REQUIRED)

add_library(Utility INTERFACE)
target_include_directories(Utility INTERFACE)

add_subdirectory(sources/Hubbard/Helper)
add_subdirectory(sources/Handler)
add_subdirectory(sources/Hubbard/Models)
add_subdirectory(sources/Hubbard/DensityOfStates)

target_link_libraries(HubbardMeanField PRIVATE 
    Models
    DensityOfStates
    Handler
    Helper
    SymbolicOperators
    Boost::iostreams
    Boost::serialization
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    ZLIB::ZLIB
)

target_link_options(HubbardMeanField PRIVATE "-Wl,--trace-symbol")