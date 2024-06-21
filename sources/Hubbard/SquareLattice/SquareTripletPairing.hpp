#pragma once
#include "HubbardCDW.hpp"

namespace Hubbard::SquareLattice {
	class SquareTripletPairing :
		public HubbardCDW
	{
	private:
		void init() override;
	protected:
		complex_prec tau_sc, theta_sc;

		inline global_floating_type theta(global_floating_type k_x, global_floating_type k_y) {
			return sin(k_x) - sin(k_y);
		};
		inline global_floating_type theta(const NumericalMomentum<2>& ks) {
			return theta(ks[0], ks[1]);
		};

		virtual void fillHamiltonian(const NumericalMomentum<2>& k_values) override;
		virtual void addToParameterSet(ParameterVector& F, const NumericalMomentum<2>& k_values) override;
	public:
		SquareTripletPairing(const ModelParameters& _params);

		virtual ModelAttributes<global_floating_type> computePhases() override;
	};
}