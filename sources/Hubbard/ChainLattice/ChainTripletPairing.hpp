#pragma once
#include "Model1D.hpp"

namespace Hubbard::ChainLattice {
	class ChainTripletPairing :
		public Model1D<complex_prec>
	{
	private:
	protected:

		virtual void fillHamiltonian(const NumericalMomentum<1>& k_x) override;

		virtual void addToParameterSet(const SpinorMatrix& rho, ParameterVector& F, const NumericalMomentum<1>& k_x) override;
	public:
		explicit ChainTripletPairing(const ModelParameters& _params);
		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy=PhaseDebuggingPolicy{}) override;
	};
}