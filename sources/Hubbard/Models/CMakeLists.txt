# CMake for Models
# requires: utility, symbolic_operators, DensityOfStates

add_library(Models 
    EMCoupling.cpp
    ModelParameters.cpp
    SquareLattice/HubbardCDW.cpp
    SquareLattice/SquareTripletPairing.cpp
    SquareLattice/UsingBroyden.cpp
    ChainLattice/ChainTripletPairing.cpp
)
SET_COMPILER_FLAGS(Models)
target_include_directories(Models PUBLIC DensityOfStates)
target_link_libraries(Models PUBLIC DensityOfStates)