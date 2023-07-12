#pragma once
#include "HubbardCDW.hpp"

namespace Hubbard::SquareLattice {
	class SquareTripletPairing :
		public HubbardCDW
	{
	protected:
		complex_prec tau_sc, theta_sc;

		inline double tau(double k_x, double k_y) {
			return sin(k_x) + sin(k_y);
		};
		inline double theta(double k_x, double k_y) {
			return sin(k_x) - sin(k_y);
		};

		virtual void fillHamiltonianHelper(va_list args) override;

		inline void addToParameterSetHelper(const SpinorMatrix& rho, ParameterVector& F, va_list args) override {
			UNPACK_2D;
			HubbardCDW::addToParameterSet(rho, F, 2, k_x, k_y);

			F(16) -= tau(k_x, k_y) * rho(6, 2);
			F(17) -= theta(k_x, k_y) * rho(6, 2);
		};
	public:
		SquareTripletPairing(const ModelParameters& _params);
		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy=PhaseDebuggingPolicy{}) override;
	};
}