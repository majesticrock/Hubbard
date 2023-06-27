#include "UsingBroyden.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "../../Utility/Roots_Broyden.hpp"

constexpr int NUMBER_OF_PARAMETERS = 8;

namespace Hubbard::SquareLattice {
	void UsingBroyden::init()
	{
		this->hamilton = SpinorMatrix::Zero(4, 4);

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
		double_prec eps = renormalizedEnergy_up(GAMMA);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = renormalizedEnergy_down(GAMMA);
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}

	UsingBroyden::UsingBroyden(const ModelParameters& _params)
		: Model(_params)
	{
		init();
	}

	UsingBroyden::UsingBroyden(const ModelParameters& _params, const BaseAttributes& startingValues)
		: Model(_params, startingValues)
	{
		init();
	}

	BaseModelRealAttributes UsingBroyden::computePhases(const bool print/*=false*/)
	{
		std::function<void(const ParameterVector&, ParameterVector&)> func = [&](const ParameterVector& x, ParameterVector& F) {
			iterationStep(x, F);
		};
		ParameterVector f0 = ParameterVector::Zero(NUMBER_OF_PARAMETERS);
		for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
		{
			f0(i) = *(parameterMapper[i]);
		}
		ParameterVector x0 = f0;

		for (size_t i = 0; i < 300 && f0.squaredNorm() > 1e-15; i++)
		{
			func(x0, f0);
			for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
			{
				x0(i) = *(parameterMapper[i]);
			}

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

		return BaseModelRealAttributes(*this);
	}
}