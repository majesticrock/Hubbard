#include "UnknownBoundaryHandler.hpp"

#ifdef _DEBUG
#define _NO_MPI
#endif
#ifndef _NO_MPI
#include <mpi.h>
#endif

#include <chrono>
#include <vector>
#include "Hubbard/Constants.hpp"

using data_vector = std::vector<Hubbard::global_floating_type>;
constexpr size_t MAX_ITERATIONS = 50;

void UnknownBoundaryHandler::execute(Utility::InputFileReader& input) const {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    data_vector local_data(FIRST_IT_STEPS);
    data_vector recieve_data(GLOBAL_IT_STEPS);
    
    for(int i = 0; i < FIRST_IT_STEPS; ++i){
        std::unique_ptr<Hubbard::Helper::ModeHelper> modeHelper{getHelper(input)};
        if(modeHelper->matrix_is_negative()){
            // Computed phase is not the thermal equilibrium
            std::cout << "Matrix is negative" << std::endl;
        } else{
            // Computed phase is the thermal equilibrium
            std::cout << "Matrix is positive" << std::endl;
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Unknown boundary computation time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
}