#pragma once
#include "Model.hpp"
#include "../BaseModelAttributes.hpp"

namespace Hubbard::SquareLattice {
	class HubbardCDW : public Model<complex_prec>
	{
	private:
		void init();
	protected:
		complex_prec gamma_cdw, xi_cdw;
		complex_prec gamma_afm, xi_afm;
		complex_prec gamma_eta, xi_eta;
		complex_prec xi_occupation_up, xi_occupation_down;

		virtual void fillHamiltonianHelper(va_list args) override;

		inline void addToParameterSetHelper(const SpinorMatrix& rho, ParameterVector& F, va_list args) override {
			UNPACK_2D;

			F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
			F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
			F(2) -= (rho(0, 2) + rho(1, 3)); // SC
			F(3) -= gamma(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Gamma SC
			F(4) -= xi(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Xi SC
			F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
			F(6) -= gamma(k_x, k_y) * (rho(0, 0) - rho(1, 1)).real(); // Gamma Occupation Up
			F(7) += gamma(k_x, k_y) * (rho(2, 2) - rho(3, 3)).real(); // Gamma Occupation Down
			F(8) += I * gamma(k_x, k_y) * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Gamma CDW
			F(9) += I * xi(k_x, k_y) * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Xi CDW
			F(10) += I * gamma(k_x, k_y) * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Gamma AFM
			F(11) += I * xi(k_x, k_y) * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Xi AFM
			F(12) -= xi(k_x, k_y) * (rho(0, 0) - rho(1, 1)).real(); // Xi Occupation Up
			F(13) += xi(k_x, k_y) * (rho(2, 2) - rho(3, 3)).real(); // Xi Occupation Down
			F(14) -= gamma(k_x, k_y) * (rho(0, 3) - rho(1, 2)); // Gamma eta
			F(15) -= xi(k_x, k_y) * (rho(0, 3) - rho(1, 2)); // Xi eta
		};
	public:
		HubbardCDW(const ModelParameters& _params);
		HubbardCDW(const ModelParameters& _params, const BaseAttributes& startingValues);
		virtual BaseModelRealAttributes computePhases(const bool print = false) override;
	};
}