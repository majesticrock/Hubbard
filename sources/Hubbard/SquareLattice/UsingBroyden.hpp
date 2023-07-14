#pragma once
#include "Model2D.hpp"

namespace Hubbard::SquareLattice {
	class UsingBroyden : public Model2D<double>
	{
	private:
		void init();
		const size_t _MaxPreBroydenIterations;
	protected:
		virtual void fillHamiltonian(const NumericalMomentum<2>& k_values) override;

		virtual void addToParameterSet(const SpinorMatrix& rho, ComplexParameterVector& F, const NumericalMomentum<2>& k_values) override;
	public:
		UsingBroyden(const ModelParameters& _params, size_t MaxPreBroydenIterations = 300U);
		UsingBroyden(const ModelParameters& _params, const BaseAttributes& startingValues, size_t MaxPreBroydenIterations = 300U);

		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy = PhaseDebuggingPolicy{}) override;
	};
}