#include "UsingBroyden.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "../Utility/Roots_Broyden.hpp"

namespace Hubbard {
	void UsingBroyden::fillHamiltonian(double_prec k_x, double_prec k_y)
	{
		hamilton.fill(0);
		const double_prec GAMMA = gamma(k_x, k_y);
		const double_prec XI = xi(k_x, k_y);

		hamilton(0, 1) = delta_cdw - delta_afm;;
		hamilton(0, 2) = delta_sc + I * (gamma_sc * GAMMA + xi_sc * XI);
		hamilton(0, 3) = I * delta_eta;

		hamilton(1, 2) = I * delta_eta;
		hamilton(1, 3) = delta_sc - I * (gamma_sc * GAMMA + xi_sc * XI);
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

	UsingBroyden::UsingBroyden(const ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: Model(_params, _number_of_basis_terms, _start_basis_at)
	{
		this->delta_cdw = (std::abs(U) + V) * 0.5 + 0.1;
		this->delta_sc = std::abs(U + std::abs(V)) * 0.3 + 0.05;
		if (V > 0) {
			this->delta_sc *= 0;
		}
		else if (V < 0) {
			this->delta_cdw *= 0;
		}
		if (U > 0) {
			this->delta_afm = std::abs(U - std::abs(V)) * 0.5 + 0.1;
		}
		else {
			this->delta_afm = 0;
		}
		this->delta_eta = 0;//U * 0.1;
		this->delta_occupation_up = V * 0.2;
		this->delta_occupation_down = V * 0.2;
		this->gamma_sc = V * 0.05;
		this->xi_sc = std::abs(V) * 0.2;

		this->hamilton = SpinorMatrix::Zero(4, 4);
	}

	Model::data_set UsingBroyden::computePhases(const bool print/*=false*/)
	{
		SpinorMatrix rho = SpinorMatrix::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

		auto lambda_func = [&](const ParameterVector& x, ParameterVector& F) {
			F.fill(0);
			delta_cdw = x(0);
			delta_afm = x(1);
			delta_sc = x(2);
			gamma_sc = x(3);
			xi_sc = x(4);
			delta_eta = x(5);
			delta_occupation_up = x(6);
			delta_occupation_down = x(7);

			complex_prec c_sc = { 0, 0 }, c_eta = { 0, 0 };
			complex_prec c_gamma_sc = { 0, 0 }, c_xi_sc = { 0, 0 };

			F = ParameterVector::Zero();
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * L_PI) / (Constants::K_DISCRETIZATION);
				for (int l = -Constants::K_DISCRETIZATION; l < 0; l++)
				{
					double_prec k_y = (l * L_PI) / (Constants::K_DISCRETIZATION);

					fillHamiltonian(k_x, k_y);
					solver.compute(hamilton);
					rho.fill(0);
					for (int i = 0; i < rho.rows(); i++)
					{
						rho(i, i) = 1 - fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * solver.eigenvectors().adjoint();

					F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
					F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
					c_sc -= (rho(0, 2) + rho(1, 3)); // SC
					c_gamma_sc -= gamma(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Gamma SC
					c_xi_sc -= xi(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Xi SC
					c_eta -= (rho(0, 3) + rho(1, 2)); // Eta
					F(6) -= cos(k_x) * (rho(0, 0) - rho(1, 1)).real(); // Occupation Up
					F(7) += cos(k_x) * (rho(2, 2) - rho(3, 3)).real(); // Occupation Down
				}
			}

			{ // Checks for numerical accurarcy
				const double ERROR_MARGIN = 1e-10 * BASIS_SIZE;
				if (std::abs(c_sc.imag()) > ERROR_MARGIN) {
					std::cout << "sc: " << c_sc << std::endl;
				}
				if (std::abs(c_eta.real()) > ERROR_MARGIN) {
					std::cout << "eta: " << c_eta << std::endl;
				}
				if (std::abs(c_gamma_sc.real()) > ERROR_MARGIN) {
					std::cout << "xi sc x: " << c_gamma_sc << std::endl;
				}
				if (std::abs(c_xi_sc.real()) > ERROR_MARGIN) {
					std::cout << "xi sc y: " << c_xi_sc << std::endl;
				}
			}

			F(2) = c_sc.real();
			F(3) = c_gamma_sc.imag();
			F(4) = c_xi_sc.imag();
			F(5) = c_eta.imag();

			F(0) *= 0.5 * U_OVER_N - 4 * V_OVER_N; // CDW
			F(1) *= 0.5 * U_OVER_N; // AFM
			F(2) *= U_OVER_N; // SC
			F(3) *= V_OVER_N; // Gamma SC
			F(4) *= V_OVER_N; // Xi SC
			F(5) *= U_OVER_N; // Eta
			F(6) *= V_OVER_N; // Occupation Up
			F(7) *= V_OVER_N; // Occupation Down

			for (size_t i = 0; i < F.size(); i++)
			{
				if (std::abs(F(i)) < 1e-14) {
					F(i) = 0;
				}
			}

			setParameters(F);
			F -= x;
		};
		std::function<void(const ParameterVector&, ParameterVector&)> func = lambda_func;
		ParameterVector f0;
		f0 << delta_cdw, delta_afm, delta_sc, gamma_sc, xi_sc, delta_eta, delta_occupation_up, delta_occupation_down;
		ParameterVector x0;
		x0 << delta_cdw, delta_afm, delta_sc, gamma_sc, xi_sc, delta_eta, delta_occupation_up, delta_occupation_down;
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

		data_set ret;
		Utility::NumericalSolver::Roots::Broyden<double_prec, 8> broyden_solver;
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