#pragma once
#include "../Selfconsistency/BroydenSolver.hpp"
#include "DOSBasedModel.hpp"

namespace Hubbard::DOSModels {
	template <class DOS>
	class BroydenDOS : public DOSBasedModel<double, DOS> {
	private:
		const size_t _MaxPreBroydenIterations;

	public:
		explicit BroydenDOS(const ModelParameters& _params, size_t MaxPreBroydenIterations = 300U)
			: DOSBasedModel<double, DOS>(_params), _MaxPreBroydenIterations(MaxPreBroydenIterations) {};
		BroydenDOS(const ModelParameters& _params, const ModelAttributes<double>& startingValues, size_t MaxPreBroydenIterations = 300U)
			: DOSBasedModel<double, DOS>(_params, startingValues), _MaxPreBroydenIterations(MaxPreBroydenIterations) {};

		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy = PhaseDebuggingPolicy{}) override
		{
			Selfconsistency::BroydenSolver solver(this, &this->model_attributes, _MaxPreBroydenIterations);
			return solver.computePhases(debugPolicy);
		};
	};
}