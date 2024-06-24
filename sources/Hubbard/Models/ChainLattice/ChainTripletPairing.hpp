#pragma once
#include "Model1D.hpp"

namespace Hubbard::Models::ChainLattice {
	class ChainTripletPairing :
		public Model1D<complex_prec>
	{
	private:
		void init() override;
	protected:

		virtual void fillHamiltonian(const NumericalMomentum<1>& k_x) override;

		virtual void addToParameterSet(ParameterVector& F, const NumericalMomentum<1>& k_x) override;
	public:
		explicit ChainTripletPairing(const ModelParameters& _params);

		virtual ModelAttributes<global_floating_type> computePhases() override;
	};
}