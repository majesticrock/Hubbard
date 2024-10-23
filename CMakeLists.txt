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

option(USE_MKL "Use Intel MKL for LAPACK" OFF)
if(USE_MKL)
    #find_package(MKL REQUIRED)
    #if(MKL_FOUND)
    #    message(STATUS "Using Intel MKL for LAPACK")
    #    set(BLAS_LIBRARIES ${MKL_LIBRARIES})
    #    set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
    #    add_definitions(-DEIGEN_USE_LAPACKE)  # Enable Eigen's LAPACK support
    #else()
    #    message(FATAL_ERROR "MKL requested but not found")
    #endif()
endif()

if(USE_MKL)
    #set(MKL_THREADING "SEQ")  # You can also choose "TBB" or "OMP" for multithreading
    #set(MKL_INTERFACE "LP64") # Or "ILP64" depending on your integer size preference
    #include_directories(${MKL_INCLUDE_DIR})

    set(MKL_LIBRARIES
        -lmkl_intel_lp64
        -lmkl_sequential
        -lmkl_core
        #${MKLROOT}/lib/intel64/libmkl_intel_${MKL_INTERFACE}.a
        #${MKLROOT}/lib/intel64/libmkl_sequential.a
        #${MKLROOT}/lib/intel64/libmkl_core.a
        #-lpthread
        #-lm
        #-ldl
    )
endif()


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