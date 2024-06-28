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
SET_COMPILER_FLAGS(ContinuumSystem)

add_subdirectory(sources/Handler)

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

find_package(Boost REQUIRED COMPONENTS serialization iostreams)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost libraries not found")
endif()
target_include_directories(HubbardMeanField PRIVATE ${Boost_INCLUDE_DIRS})

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

target_link_libraries(HubbardMeanField PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX ${Boost_LIBRARIES} ZLIB::ZLIB Handler Helper Models)

