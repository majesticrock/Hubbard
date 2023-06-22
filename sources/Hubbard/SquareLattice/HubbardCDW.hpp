#pragma once
#include "Model.hpp"
#include "../BaseModelAttributes.hpp"

namespace Hubbard::SquareLattice {
	class HubbardCDW : public Model<complex_prec>, public BaseModelComplexAttributes
	{
	protected:
		complex_prec gamma_cdw, xi_cdw;
		complex_prec gamma_afm, xi_afm;
		complex_prec gamma_eta, xi_eta;

		complex_prec delta_occupation_up_y, delta_occupation_down_y;

		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const {
			return -2. * (
				(delta_occupation_up + 1.) * cos(k_x) + (delta_occupation_up_y + 1.) * cos(k_y)
				).real();
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const {
			return -2. * (
				(delta_occupation_down + 1.) * cos(k_x) + (delta_occupation_down_y + 1.) * cos(k_y)
				).real();
		};

		virtual void fillHamiltonianHelper(va_list args) override;

		inline void addToParameterSetHelper(const SpinorMatrix& rho, ParameterVector& F, va_list args) override {
			UNPACK_2D;

			F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
			F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
			F(2) -= (rho(0, 2) + rho(1, 3)); // SC
			F(3) -= gamma(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Gamma SC
			F(4) -= xi(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Xi SC
			F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
			F(6) -= cos(k_x) * (rho(0, 0) - rho(1, 1)).real(); // Occupation Up
			F(7) += cos(k_x) * (rho(2, 2) - rho(3, 3)).real(); // Occupation Down
			F(8) += I * gamma(k_x, k_y) * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Gamma CDW
			F(9) += I * xi(k_x, k_y) * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Xi CDW
			F(10) += I * gamma(k_x, k_y) * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Gamma AFM
			F(11) += I * xi(k_x, k_y) * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Xi AFM
			F(12) -= cos(k_y) * (rho(0, 0) - rho(1, 1)).real(); // Occupation Up y
			F(13) += cos(k_y) * (rho(2, 2) - rho(3, 3)).real(); // Occupation Down y
			F(14) -= gamma(k_x, k_y) * (rho(0, 3) - rho(1, 2)); // Gamma eta
			F(15) -= xi(k_x, k_y) * (rho(0, 3) - rho(1, 2)); // Xi eta
		};
	public:
		HubbardCDW(const ModelParameters& _params);
		virtual PhaseDataSet computePhases(const bool print = false) override;
	};
}