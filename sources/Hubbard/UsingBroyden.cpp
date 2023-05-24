#include "UsingBroyden.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "../Utility/Roots_Broyden.hpp"

namespace Hubbard {
	void UsingBroyden::computeChemicalPotential()
	{
		Model::computeChemicalPotential();
		chemical_potential += 4 * V;
	}
	void UsingBroyden::fillHamiltonian(double_prec k_x, double_prec k_y)
	{
		hamilton.fill(0);

		hamilton(0, 1) = delta_cdw - delta_afm;;
		hamilton(0, 2) = delta_sc + I * (gamma_sc * gamma(k_x, k_y) + xi_sc * xi(k_x, k_y));
		hamilton(0, 3) = I * delta_eta;

		hamilton(1, 2) = I * delta_eta;
		hamilton(1, 3) = delta_sc - I * (gamma_sc * gamma(k_x, k_y) + xi_sc * xi(k_x, k_y));
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

	UsingBroyden::UsingBroyden(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: Model(_params, _number_of_basis_terms, _start_basis_at), V(_params.V)
	{
		this->delta_cdw = std::abs(U - V) * 0.5 + 0.1;
		this->delta_sc = std::abs(U + V) * 0.5 + 0.1;
		if (V > 0) {
			this->delta_sc *= 0.25;
		}
		else if (V < 0) {
			this->delta_cdw *= 0;
		}
		this->delta_afm = (U > 0 && U > 4*V) ? -this->delta_cdw : this->delta_cdw;
		
		this->delta_eta = 0;//std::abs(U + V) * 0.2;
		this->delta_occupation_up = V * 0.1;
		this->delta_occupation_down = std::abs(V) * 0.1;
		this->gamma_sc = std::abs(V) * 0.5;
		this->xi_sc = std::abs(V) * 0.5;

		this->V_OVER_N = V / BASIS_SIZE;

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
					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * (solver.eigenvectors().adjoint());

					F(0) += (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real();
					F(1) += (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real();
					c_sc       += (rho(2, 0) + rho(3, 1));
					c_gamma_sc += gamma(k_x, k_y) * (rho(2, 0) - rho(3, 1));
					c_xi_sc    += xi(k_x, k_y)    * (rho(2, 0) - rho(3, 1));
					c_eta      += (rho(0, 3) + rho(1, 2));
					F(6)       += cos(k_x) * (rho(0, 0) - rho(1, 1)).real();
					F(7)       -= cos(k_x) * (rho(2, 2) - rho(3, 3)).real();
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

			setParameters(F);
			F -= x;
		};
		std::function<void(const ParameterVector&, ParameterVector&)> func = lambda_func;
		ParameterVector f0;
		f0 << delta_cdw, delta_afm, delta_sc, gamma_sc, xi_sc, delta_eta, delta_occupation_up, delta_occupation_down;
		ParameterVector x0;
		x0 << delta_cdw, delta_afm, delta_sc, gamma_sc, xi_sc, delta_eta, delta_occupation_up, delta_occupation_down;
		for (size_t i = 0; i < 200 && f0.squaredNorm() > 1e-15; i++)
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

			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow(x0);
			}
		}

		Utility::NumericalSolver::Roots::Broyden<double_prec, 8> broyden_solver;
		if (!broyden_solver.compute(func, x0)) {
			std::cerr << "No convergence for [T, U, V] = [" << this->temperature << ", " << this->U << ", " << this->V << "]" << std::endl;
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

		data_set ret;
		ret.delta_cdw = delta_cdw;
		ret.delta_afm = delta_afm;
		ret.delta_sc = delta_sc;
		ret.gamma_sc = gamma_sc;
		ret.xi_sc = xi_sc;
		ret.delta_eta = delta_eta;

		return ret;
	}
}