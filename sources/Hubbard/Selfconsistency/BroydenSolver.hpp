#pragma once
#include "IterativeSolver.hpp"

namespace Hubbard::Selfconsistency {
	class BroydenSolver : public IterativeSolver<double>
	{
	private:
		const size_t _MaxPreBroydenIterations;
	public:
		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy& debugPolicy);

		BroydenSolver() = delete;
		BroydenSolver(BaseModel<double>* model_ptr, ModelAttributes<double>* attribute_ptr, size_t MaxPreBroydenIterations)
			: IterativeSolver<double>(model_ptr, attribute_ptr), _MaxPreBroydenIterations(MaxPreBroydenIterations) {};
	};
}