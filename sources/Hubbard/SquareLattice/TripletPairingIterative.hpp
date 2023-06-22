#pragma once
#include "HubbardCDW.hpp"

namespace Hubbard::SquareLattice {
	class TripletPairingIterative :
		public HubbardCDW
	{
	protected:
		complex_prec tau_sc, theta_sc;

		inline double_prec tau(double_prec k_x, double_prec k_y) {
			return sin(k_x) + sin(k_y);
		};
		inline double_prec theta(double_prec k_x, double_prec k_y) {
			return sin(k_x) - sin(k_y);
		};

		virtual void fillHamiltonianHelper(va_list args) override;

		inline void addToParameterSetHelper(const SpinorMatrix& rho, ParameterVector& F, va_list args) override {
			UNPACK_2D;
			HubbardCDW::addToParameterSet(rho, F, k_x, k_y);

			F(16) -= tau(k_x, k_y) * rho(6, 2);
			F(17) -= theta(k_x, k_y) * rho(6, 2);
		};
	public:
		TripletPairingIterative(const ModelParameters& _params);
		virtual PhaseDataSet computePhases(const bool print = false) override;
	};
}