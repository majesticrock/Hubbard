#pragma once
#include "Model.hpp"
#include "../BaseModelAttributes.hpp"

namespace Hubbard::SquareLattice {
	class UsingBroyden : public Model<double_prec>, public BaseModelRealAttributes
	{
	protected:
		virtual void fillHamiltonianHelper(va_list args) override;

		inline void addToParameterSetHelper(const SpinorMatrix& rho, ComplexParameterVector& F, va_list args) override {
			UNPACK_2D;

			F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
			F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
			F(2) -= (rho(0, 2) + rho(1, 3)); // SC
			F(3) -= gamma(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Gamma SC
			F(4) -= xi(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Xi SC
			F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
			F(6) -= cos(k_x) * (rho(0, 0) - rho(1, 1)).real(); // Occupation Up
			F(7) += cos(k_x) * (rho(2, 2) - rho(3, 3)).real(); // Occupation Down
		};

		virtual inline void complexParametersToReal(const ComplexParameterVector& c, ParameterVector& r) const override {
			{ // Checks for numerical accurarcy
				const double ERROR_MARGIN = 1e-10 * Constants::BASIS_SIZE;
				if (std::abs(c(2).imag()) > ERROR_MARGIN) {
					std::cout << "sc: " << c(2) << std::endl;
				}
				if (std::abs(c(3).imag()) > ERROR_MARGIN) {
					std::cout << "gamma sc: " << c(3) << std::endl;
				}
				if (std::abs(c(3).real()) > ERROR_MARGIN) {
					std::cout << "xi sc y: " << c(4) << std::endl;
				}
				if (std::abs(c(5).real()) > ERROR_MARGIN) {
					std::cout << "eta: " << c(5) << std::endl;
				}
			}

			r(0) = c(0).real(); // CDW
			r(1) = c(1).real(); // AFM
			r(2) = c(2).real(); // SC
			r(3) = c(3).real(); // Gamma SC
			r(4) = c(4).imag(); // Xi SC
			r(5) = c(5).imag(); // Eta
			r(6) = c(6).real(); // Occupation Up
			r(7) = c(7).real(); // Occupation Down
		};
	public:
		UsingBroyden(const ModelParameters& _params);

		PhaseDataSet computePhases(const bool print = false) override;
	};
}