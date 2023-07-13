#pragma once
#include "BaseModel.hpp"
#include <memory>

namespace Hubbard::Selfconsistency{
    class IterativeSolver{
	private:
		BaseModel<double>* real_model{};
		BaseModel<std::complex<double>>* complex_model{};
		const bool isReal{};
	public:
		ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy=PhaseDebuggingPolicy{});

		IterativeSolver(BaseModel<double>* model_ptr);
		IterativeSolver(BaseModel<std::complex<double>>* model_ptr);
    };
}