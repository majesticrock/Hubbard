# CMake for DensityOfStates
# requires: Boost, openMP, Utility

add_library(DensityOfStates 
    BaseDOS.cpp
    SimpleCubic.cpp
    Square.cpp
)
SET_COMPILER_FLAGS(DensityOfStates)
target_include_directories(DensityOfStates PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(DensityOfStates PUBLIC ZLIB::ZLIB ${Boost_LIBRARIES})