#pragma once
#include "Model2D.hpp"

namespace Hubbard::SquareLattice {
	class UsingBroyden : public Model2D<double>
	{
	private:
		const size_t MaxPreBroydenIterations;
		void init();
	protected:
		virtual void fillHamiltonian(const std::array<double, 2>& k_values) override;

		virtual void addToParameterSetHelper(const SpinorMatrix& rho, ComplexParameterVector& F, const std::array<double, 2>& k_values) override;
	public:
		explicit UsingBroyden(const ModelParameters& _params, size_t _MaxPreBroydenIterations = 300U);
		UsingBroyden(const ModelParameters& _params, const BaseAttributes& startingValues, size_t _MaxPreBroydenIterations = 300U);

		ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy=PhaseDebuggingPolicy{}) override;
	};
}