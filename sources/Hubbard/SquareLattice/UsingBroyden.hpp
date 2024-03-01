#pragma once
#include "Model2D.hpp"

namespace Hubbard::SquareLattice {
	class UsingBroyden : public Model2D<global_floating_type>
	{
	private:
		void init() override;
		const size_t _MaxPreBroydenIterations;
	protected:
		virtual void fillHamiltonian(const NumericalMomentum<2>& k_values) override;

		virtual void addToParameterSet(ComplexParameterVector& F, const NumericalMomentum<2>& k_values) override;
	public:
		UsingBroyden(const ModelParameters& _params, size_t MaxPreBroydenIterations = 300U);
		UsingBroyden(const ModelParameters& _params, const BaseAttributes& startingValues, size_t MaxPreBroydenIterations = 300U);

		virtual ModelAttributes<global_floating_type> computePhases(const PhaseDebuggingPolicy debugPolicy = WarnNoConvergence) override;
	};
}