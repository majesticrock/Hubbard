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
    target_link_libraries(Helper PUBLIC ${MKL_LIBRARIES})
endif()
target_include_directories(Helper PUBLIC ${SymbolicOperators_INCLUDE_DIRS})
target_link_libraries(Helper PUBLIC OpenMP::OpenMP_CXX SymbolicOperators::SymbolicOperators Models)