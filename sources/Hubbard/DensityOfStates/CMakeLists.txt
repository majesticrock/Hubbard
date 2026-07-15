# CMake for DensityOfStates
# requires: Boost, openMP, utility

add_library(DensityOfStates 
    BaseDOS.cpp
    SimpleCubic.cpp
    Square.cpp
)

target_link_libraries(DensityOfStates 
    PUBLIC 
        mrock::utility 
        OpenMP::OpenMP_CXX
        hubbard_options
)