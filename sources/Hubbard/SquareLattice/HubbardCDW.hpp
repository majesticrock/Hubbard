#pragma once
#include "Model2D.hpp"

namespace Hubbard::SquareLattice {
	class HubbardCDW : public Model2D<complex_prec>
	{
	private:
		void init();
	protected:
		virtual void fillHamiltonian(const NumericalMomentum<2>& k_values) override;

		virtual void addToParameterSet(ParameterVector& F, const NumericalMomentum<2>& k_values) override;
	public:
		HubbardCDW(const ModelParameters& _params);

		template<typename StartingValuesDataType>
		HubbardCDW(const ModelParameters& _params, const ModelAttributes< StartingValuesDataType >& startingValues)
			: Model2D(_params, startingValues) {};

		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy = PhaseDebuggingPolicy{}) override;
	};
}