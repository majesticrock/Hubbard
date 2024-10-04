set(HELPER_SOURCES
    GeneralBasis.cpp
    ModeHelper.cpp
    PhaseHelper.cpp
    Plaquette.cpp
    SquareGeneral.cpp
    SquareXP.cpp
    TermOnSquare.cpp
    XPModes.cpp
)

#add_subdirectory(../Models ${CMAKE_BINARY_DIR}/Models)
#add_subdirectory(../DensityOfStates ${CMAKE_BINARY_DIR}/DensistyOfStates)

add_library(Helper STATIC ${HELPER_SOURCES})
SET_COMPILER_FLAGS(Helper)
target_include_directories(Helper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#target_link_libraries(Helper PUBLIC Models PUBLIC DensityOfStates)
