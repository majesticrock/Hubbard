#pragma once
#include "IterativeSolver.hpp"

namespace Hubbard::Selfconsistency {
	class BroydenSolver : public IterativeSolver<global_floating_type>
	{
	private:
		const size_t _MaxPreBroydenIterations;
	public:
		virtual const ModelAttributes<global_floating_type>& computePhases(const PhaseDebuggingPolicy& debugPolicy);

		BroydenSolver() = delete;
		BroydenSolver(BaseModel<global_floating_type>* model_ptr, ModelAttributes<global_floating_type>* attribute_ptr, size_t MaxPreBroydenIterations)
			: IterativeSolver<global_floating_type>(model_ptr, attribute_ptr), _MaxPreBroydenIterations(MaxPreBroydenIterations) {};
	};
}