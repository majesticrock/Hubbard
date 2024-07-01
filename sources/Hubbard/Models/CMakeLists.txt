set(MODEL_SOURCES
    EMCoupling.cpp
    ModelParameters.cpp
    SquareLattice/HubbardCDW.cpp
    SquareLattice/SquareTripletPairing.cpp
    SquareLattice/UsingBroyden.cpp
    ChainLattice/ChainTripletPairing.cpp
)

add_library(Models STATIC ${MODEL_SOURCES})
target_include_directories(Models PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
SET_COMPILER_FLAGS(Models)

find_package(SymbolicOperators REQUIRED CONFIG)
if(NOT SymbolicOperators_FOUND)
    message(FATAL_ERROR "SymbolicOperators library not found")
endif()
target_include_directories(Models PUBLIC ${SymbolicOperators_INCLUDE_DIRS})
target_link_libraries(Models PUBLIC SymbolicOperators::SymbolicOperators)

find_package(Utility REQUIRED)
if(NOT Utility_FOUND)
    message(Utility_ERROR "Utility library not found")
endif()
target_include_directories(Models PUBLIC ${Utility_INCLUDE_DIRS})
target_link_libraries(Models PUBLIC Utility::Utility)