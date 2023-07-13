#include "IterativeSolver.hpp"

namespace Hubbard::Selfconsistency
{
    IterativeSolver::IterativeSolver(BaseModel<double>* model_ptr)
			: real_model{model_ptr}, isReal{true} {};

    IterativeSolver::IterativeSolver(BaseModel<std::complex<double>>* model_ptr)
			: complex_model{model_ptr}, isReal{false} {};

    ModelAttributes<double> IterativeSolver::computePhases(const PhaseDebuggingPolicy debugPolicy/*=PhaseDebuggingPolicy{}*/)
	{
		constexpr double EPSILON = 1e-12;
		double error = 100;
		constexpr size_t MAX_STEPS = 1000;
        const size_t NUMBER_OF_PARAMETERS = model_attributes.size();

       // decltype(auto) model = (isReal ? real_model : complex_model);

		ParameterVector f0{ ParameterVector::Zero(NUMBER_OF_PARAMETERS) };
		std::copy(model->model_attributes.selfconsistency_values.begin(), model->model_attributes.selfconsistency_values.end(), f0.begin());

		ParameterVector x0 = f0;

		if (debugPolicy.printAll) {
			std::cout << "-1:\t" << std::fixed << std::setprecision(8);
			printAsRow<-1>(x0);
		}
		model->model_attributes.converged = true;
		for (size_t i = 0U; i < MAX_STEPS && error > EPSILON; ++i)
		{
			model->iterationStep(x0, f0);

            for (size_t j = 0U; j < NUMBER_OF_PARAMETERS; ++j)
			{
				if(std::abs(x0[j]) > 1e-10){
					if(std::abs((x0[j] + model_attributes[j]) / x0[j]) < 1e-12){
						// Sign flipping behaviour
						if (debugPolicy.convergenceWarning){
							std::cerr << "Sign flipper for [T U V] = [" << std::fixed << std::setprecision(8)
							<< this->temperature << " " << this->U << " " << this->V << "]" << std::endl;
						}

						std::fill(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), 0.);
						model_attributes.converged = false;
						return ModelAttributes<double>(this->model_attributes);
					}
				}
			}

			error = f0.norm();
			std::copy(model->model_attributes.selfconsistency_values.begin(), model->model_attributes.selfconsistency_values.end(), x0.begin());

			if (debugPolicy.printAll) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow<-1>(x0);
			}
			if (i == MAX_STEPS - 1) {
				if (debugPolicy.convergenceWarning){
					std::cerr << "No convergence for [T U V] = [" << std::fixed << std::setprecision(8)
					<< this->temperature << " " << this->U << " " << this->V << "]" << std::endl;
				}

				std::fill(model->model_attributes.selfconsistency_values.begin(), model->model_attributes.selfconsistency_values.end(), 0.);
				model->model_attributes.converged = false;
			}
		}

		return ModelAttributes<double>(model->model_attributes);
	}
}