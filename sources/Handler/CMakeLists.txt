# CMake for Handler
# requires: MPI, Utility, Models, DensityOfStates, Helper

add_library(Handler 
    HandlerBase.cpp
    ModeDispersionHandler.cpp
    ModeHandler.cpp
    PhaseHandler.cpp
    TestHandler.cpp
    UnknownBoundaryHandler.cpp
)
SET_COMPILER_FLAGS(Handler)
target_include_directories(Handler PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(Handler PUBLIC ZLIB::ZLIB ${Boost_LIBRARIES} MPI::MPI_CXX)