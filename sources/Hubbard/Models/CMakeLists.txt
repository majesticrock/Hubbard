# CMake for Models
# requires: utility, symbolic_operators, DensityOfStates

add_library(Models 
    ModelParameters.cpp
    SquareLattice/HubbardCDW.cpp
    SquareLattice/SquareTripletPairing.cpp
    SquareLattice/UsingBroyden.cpp
    ChainLattice/ChainTripletPairing.cpp
)

target_include_directories(Models PUBLIC DensityOfStates)
target_link_libraries(Models PUBLIC DensityOfStates hubbard_options)