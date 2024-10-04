set(DOS_SOURCES
    BaseDOS.cpp
    SimpleCubic.cpp
    Square.cpp
)

add_library(DensityOfStates STATIC ${DOS_SOURCES})
target_include_directories(DensityOfStates PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
SET_COMPILER_FLAGS(DensityOfStates)

#find_package(Utility REQUIRED)
#if(NOT Utility_FOUND)
#    message(Utility_ERROR "Utility library not found")
#endif()
#target_include_directories(DensityOfStates PRIVATE ${Utility_INCLUDE_DIRS})
#target_link_libraries(DensityOfStates PUBLIC Utility::Utility)
