# CMake for Models
# requires: Utility, SymbolicOperators, DensityOfStates

add_library(Models 
    EMCoupling.cpp
    ModelParameters.cpp
    SquareLattice/HubbardCDW.cpp
    SquareLattice/SquareTripletPairing.cpp
    SquareLattice/UsingBroyden.cpp
    ChainLattice/ChainTripletPairing.cpp
)
SET_COMPILER_FLAGS(Models)
target_include_directories(Models PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(Models PUBLIC ZLIB::ZLIB ${Boost_LIBRARIES})