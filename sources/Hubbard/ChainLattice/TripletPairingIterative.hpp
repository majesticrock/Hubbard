#pragma once
#include "Model1D.hpp"

namespace Hubbard::ChainLattice {
	class TripletPairingIterative :
		public Model1D<complex_prec>
	{
	private:
	protected:
		inline double_prec tau(double_prec k_x) {
			return sin(k_x);
		};
		virtual void fillHamiltonianHelper(va_list args) override;

		inline void addToParameterSetHelper(const SpinorMatrix& rho, ParameterVector& F, va_list args) override {
			UNPACK_1D;

			F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
			F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
			F(2) -= (rho(0, 2) + rho(1, 3)); // SC
			F(3) -= gamma(k_x) * (rho(0, 2) - rho(1, 3)); // Gamma SC
			F(4) -= tau(k_x) * rho(6, 2); // Tau SC
			F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
			F(6) -= gamma(k_x) * (rho(0, 0) - rho(1, 1)).real(); // Gamma Occupation Up
			F(7) += gamma(k_x) * (rho(2, 2) - rho(3, 3)).real(); // Gamma Occupation Down
		};
	public:
		TripletPairingIterative(const ModelParameters& _params);
		virtual BaseModelRealAttributes computePhases(const bool print = false) override;
	};
}