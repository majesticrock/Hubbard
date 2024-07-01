set(HANDLER_SOURCES
    HandlerBase.cpp
    ModeDispersionHandler.cpp
    ModeHandler.cpp
    PhaseHandler.cpp
    TestHandler.cpp
    UnknownBoundaryHandler.cpp
)

add_subdirectory(../Hubbard/Helper ${CMAKE_BINARY_DIR}/Hubbard/Helper)

add_library(Handler STATIC ${HANDLER_SOURCES})
target_include_directories(Handler PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
SET_COMPILER_FLAGS(Handler)

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
target_link_libraries(Handler PUBLIC MPI::MPI_CXX Models)