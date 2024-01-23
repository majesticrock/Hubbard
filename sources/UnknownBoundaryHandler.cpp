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
constexpr size_t MAX_ITERATIONS = 30;
constexpr double ERROR_MARGIN = 1e-6;

void UnknownBoundaryHandler::execute(Utility::InputFileReader& input) const {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	double SECOND_IT_MIN = 0, SECOND_IT_MAX = input.getDouble("second_iterator_upper_limit");

	for (int i = 0; i < Hubbard::Constants::option_list.size(); i++)
	{
		if (input.getString("second_iterator_type") == Hubbard::Constants::option_list[i]) {
			SECOND_IT_MIN = model_params[i];
		}
	}
    Hubbard::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
			(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, 0.0,
			input.getString("global_iterator_type"), input.getString("second_iterator_type"));
    for (int i = 0; i < FIRST_IT_STEPS; ++i){
        modelParameters.printParameters();
        modelParameters.incrementGlobalIterator();
    }
    return;

    data_vector local_data(FIRST_IT_STEPS);
    data_vector recieve_data(GLOBAL_IT_STEPS);

    for(int i = 0; i < FIRST_IT_STEPS; ++i){
        std::unique_ptr<Hubbard::Helper::ModeHelper> modeHelper;
        double a = SECOND_IT_MIN;
        double b = SECOND_IT_MAX;
        double c = 0.5 * (a + b);

        modelParameters.setSecondIteratorExact(a);
        modeHelper = getHelper(input);
        if(modeHelper->matrix_is_negative()){
            // Computed phase is not the thermal equilibrium
            
        } else{
            // Computed phase is the thermal equilibrium
            std::cout << "Matrix is positive" << std::endl;
        }

        for(size_t iter = 0U; iter < MAX_ITERATIONS && abs(b - a) > ERROR_MARGIN; ++iter){
            modeHelper = getHelper(input);
            if(modeHelper->matrix_is_negative()){
                // Computed phase is not the thermal equilibrium
                std::cout << "Matrix is negative" << std::endl;
            } else{
                // Computed phase is the thermal equilibrium
                std::cout << "Matrix is positive" << std::endl;
            }
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Unknown boundary computation time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
}