# CMake for Handler
# requires: MPI, utility, Models, DensityOfStates, Helper

add_library(Handler 
    HandlerBase.cpp
    ModeDispersionHandler.cpp
    ModeHandler.cpp
    PhaseHandler.cpp
    TestHandler.cpp
    UnknownBoundaryHandler.cpp
)

target_link_libraries(Handler 
    PRIVATE 
        MPI::MPI_CXX 
        Models 
        Helper
        hubbard_options
)