#pragma once
#include <Utility/Selfconsistency/BroydenSolver.hpp>
#include "DOSBasedModel.hpp"

namespace Hubbard::DOSModels {
	template <class DOS>
	class BroydenDOS : public DOSBasedModel<global_floating_type, DOS> {
	private:
		const size_t _MaxPreBroydenIterations;

	public:
		explicit BroydenDOS(const ModelParameters& _params, size_t MaxPreBroydenIterations = 300U)
			: DOSBasedModel<global_floating_type, DOS>(_params), _MaxPreBroydenIterations(MaxPreBroydenIterations) {};
		BroydenDOS(const ModelParameters& _params, const ModelAttributes<global_floating_type>& startingValues, size_t MaxPreBroydenIterations = 300U)
			: DOSBasedModel<global_floating_type, DOS>(_params, startingValues), _MaxPreBroydenIterations(MaxPreBroydenIterations) {};

		virtual ModelAttributes<global_floating_type> computePhases(const PhaseDebuggingPolicy debugPolicy = PrintSteps) override
		{
			Utility::Selfconsistency::BroydenSolver<global_floating_type, BroydenDOS<DOS>, decltype(this->model_attributes)> 
				solver(this, &this->model_attributes, _MaxPreBroydenIterations);
			return solver.compute(true);
		};
	};
}