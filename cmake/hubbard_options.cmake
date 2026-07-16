# Custom architecture target
set(TARGET_ARCH
    "native"
    CACHE STRING
    "Architecture passed to compiler as -march=<arch> (empty disables)"
)

# Joined compile options
add_library(hubbard_options INTERFACE)
target_compile_features(hubbard_options
    INTERFACE
        cxx_std_20
)

target_compile_options(hubbard_options INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wextra>
)
if(TARGET_ARCH)
    target_compile_options(hubbard_options INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-march=${TARGET_ARCH}>
    )
    message(STATUS "Building for architecture ${TARGET_ARCH}!")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(hubbard_options INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-ffast-math> 
    )
else()
    target_compile_definitions(hubbard_options INTERFACE DEBUG)
endif()

#
# Use MKL if available
#
set(MKL_LINK static)
set(MKL_THREADING gnu_thread)
set(MKL_INTERFACE lp64)

find_package(MKL CONFIG QUIET HINTS $ENV{MKLROOT})

if (MKL_FOUND)
    message(STATUS "Configuring hubbard to use Intel MKL")

    target_compile_options(hubbard_options INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(hubbard_options INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(hubbard_options INTERFACE $<LINK_ONLY:MKL::MKL>)
    # Let Eigen use MKL
    target_compile_definitions(hubbard_options INTERFACE EIGEN_USE_MKL_ALL MROCK_IEOM_DO_NOT_PARALLELIZE)
else()
    message(STATUS "MKL not found or not enabled; hubbard will not use Intel MKL")
endif()