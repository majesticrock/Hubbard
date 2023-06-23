#include "UsingBroyden.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "../../Utility/Roots_Broyden.hpp"

constexpr int NUMBER_OF_PARAMETERS = 8;

namespace Hubbard::SquareLattice {
	void UsingBroyden::fillHamiltonianHelper(va_list args)
	{
		UNPACK_2D;
		hamilton.fill(0);
		const double_prec GAMMA = gamma(k_x, k_y);
		const double_prec XI = xi(k_x, k_y);

		hamilton(0, 1) = delta_cdw - delta_afm;;
		hamilton(0, 2) = delta_sc + (gamma_sc * GAMMA + I * xi_sc * XI);
		hamilton(0, 3) = I * delta_eta;

		hamilton(1, 2) = I * delta_eta;
		hamilton(1, 3) = delta_sc - (gamma_sc * GAMMA + I * xi_sc * XI);
		hamilton(2, 3) = -delta_cdw - delta_afm;

		SpinorMatrix buffer = hamilton.adjoint();
		hamilton += buffer;
		double_prec eps = renormalizedEnergy_up(k_x, k_y);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = renormalizedEnergy_down(k_x, k_y);
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}

	UsingBroyden::UsingBroyden(const ModelParameters& _params)
		: Model(_params), BaseModelRealAttributes(_params)
	{
		this->hamilton = SpinorMatrix::Zero(4, 4);

		parameterMapper = {
			&(this->delta_cdw),
			&(this->delta_afm),
			&(this->delta_sc),
			&(this->gamma_sc),
			&(this->xi_sc),
			&(this->delta_eta),
			&(this->delta_occupation_up),
			&(this->delta_occupation_down),
		};

		parameterCoefficients = {
			0.5 * U_OVER_N - 4. * V_OVER_N, // CDW
			0.5 * U_OVER_N, // AFM
			U_OVER_N, // SC
			V_OVER_N, // Gamma SC
			V_OVER_N, // Xi SC
			U_OVER_N, // Eta
			V_OVER_N, // Occupation Up
			V_OVER_N, // Occupation Down
		};
	}

	PhaseDataSet UsingBroyden::computePhases(const bool print/*=false*/)
	{
		SpinorMatrix rho = SpinorMatrix::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

		std::function<void(const ParameterVector&, ParameterVector&)> func = [&](const ParameterVector& x, ParameterVector& F) {
			iterationStep(x, F);
		};
		ParameterVector f0 = ParameterVector::Zero(NUMBER_OF_PARAMETERS);
		f0 << delta_cdw, delta_afm, delta_sc, gamma_sc, xi_sc, delta_eta, delta_occupation_up, delta_occupation_down;
		ParameterVector x0;
		x0 = f0;
		for (size_t i = 0; i < 300 && f0.squaredNorm() > 1e-15; i++)
		{
			func(x0, f0);
			x0(0) = delta_cdw;
			x0(1) = delta_afm;
			x0(2) = delta_sc;
			x0(3) = gamma_sc;
			x0(4) = xi_sc;
			x0(5) = delta_eta;
			x0(6) = delta_occupation_up;
			x0(7) = delta_occupation_down;

			for (size_t i = 0; i < x0.size(); i++)
			{
				if (std::abs(x0(i)) < 1e-14) {
					x0(i) = 0;
				}
			}

			if (print) {
				std::cout << i << ":  " << std::scientific << std::setprecision(4);
				printAsRow(x0);
			}
		}

		PhaseDataSet ret;
		Utility::NumericalSolver::Roots::Broyden<double_prec, -1> broyden_solver;
		if (!broyden_solver.compute(func, x0, 400)) {
			std::cerr << "No convergence for [T U V] = [" << std::fixed << std::setprecision(8)
				<< this->temperature << " " << this->U << " " << this->V << "]" << std::endl;
			delta_cdw = 0;
			delta_afm = 0;
			delta_sc = 0;
			gamma_sc = 0;
			xi_sc = 0;
			delta_eta = 0;
			ret.converged = false;
		}

		if (print) {
			func(x0, f0);
			std::cout << "T=" << temperature << "   U=" << U << "   V=" << V << "\n";
			std::cout << "x0 = (";
			for (int i = 0; i < x0.size(); i++)
			{
				std::cout << " " << x0(i) << " ";
			}
			std::cout << ")\nf0 = (";
			for (int i = 0; i < f0.size(); i++)
			{
				std::cout << " " << f0(i) << " ";
			}
			std::cout << ")\n -> |f0| = " << std::scientific << std::setprecision(8) << f0.norm() << std::endl;
		}

		ret.delta_cdw = delta_cdw;
		ret.delta_afm = delta_afm;
		ret.delta_sc = delta_sc;
		ret.gamma_sc = gamma_sc;
		ret.xi_sc = xi_sc;
		ret.delta_eta = delta_eta;

		return ret;
	}
}