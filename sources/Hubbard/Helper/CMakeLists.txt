# CMake for Helper
# requires: Utility, Models, SymbolicOperators, DensityOfStates

add_library(Helper 
    GeneralBasis.cpp
    ModeHelper.cpp
    PhaseHelper.cpp
    Plaquette.cpp
    SquareGeneral.cpp
    SquareXP.cpp
    TermOnSquare.cpp
    XPModes.cpp
)
SET_COMPILER_FLAGS(Helper)

if(USE_MKL)
    target_compile_options(Helper PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(Helper PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(Helper PUBLIC $<LINK_ONLY:MKL::MKL>)

    get_target_property(MKL_LINK_LIBS MKL::MKL INTERFACE_LINK_LIBRARIES)
    message(STATUS "MKL::MKL INTERFACE_LINK_LIBRARIES: ${MKL_LINK_LIBS}")

    target_compile_definitions(Helper PUBLIC EIGEN_USE_LAPACKE)
endif()

target_include_directories(Helper PUBLIC ${SymbolicOperators_INCLUDE_DIRS})
target_link_libraries(Helper PUBLIC OpenMP::OpenMP_CXX SymbolicOperators::SymbolicOperators Models)