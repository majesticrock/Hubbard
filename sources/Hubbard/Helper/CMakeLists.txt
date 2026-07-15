# CMake for Helper
# requires: utility, Models, symbolic_operators, DensityOfStates

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

target_include_directories(Helper PUBLIC ${mrock_INCLUDE_DIRS})
target_link_libraries(Helper 
    PUBLIC 
        OpenMP::OpenMP_CXX 
        mrock::symbolic_operators 
        Models
        hubbard_options
)