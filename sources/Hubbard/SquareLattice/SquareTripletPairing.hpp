#pragma once
#include "HubbardCDW.hpp"

namespace Hubbard::SquareLattice {
	class SquareTripletPairing :
		public HubbardCDW
	{
	protected:
		complex_prec tau_sc, theta_sc;

		inline double theta(double k_x, double k_y) {
			return sin(k_x) - sin(k_y);
		};
		inline double theta(const std::array<double, 2>& ks){
			return theta(ks[0], ks[1]);
		};

		virtual void fillHamiltonian(const std::array<double, 2>& k_values) override;
		virtual void addToParameterSet(const SpinorMatrix& rho, ParameterVector& F, const std::array<double, 2>& k_values) override;
	public:
		SquareTripletPairing(const ModelParameters& _params);
		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy=PhaseDebuggingPolicy{}) override;
	};
}