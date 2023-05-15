#include "UsingBroyden.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include "../Utility/Roots_Broyden.hpp"

namespace Hubbard {
	void UsingBroyden::computeChemicalPotential()
	{
		Model::computeChemicalPotential();
		chemical_potential += 0.5 * V;
	}
	void UsingBroyden::fillHamiltonian(double_prec k_x, double_prec k_y)
	{
		hamilton.fill(0);

		hamilton(0, 1) = delta_cdw_up;
		hamilton(0, 2) = delta_sc;
		hamilton(0, 3) = delta_eta;

		hamilton(1, 2) = delta_eta;
		hamilton(1, 3) = delta_sc;
		hamilton(2, 3) = -delta_cdw_down;

		Matrix_L buffer = hamilton.transpose();
		hamilton += buffer;
		const double_prec eps = renormalizedEnergy(k_x, k_y); //unperturbed_energy(k_x, k_y) - (2 * delta_occupation * (cos(k_x) + cos(k_y))); //
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}

	UsingBroyden::UsingBroyden(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: Model(_params, _number_of_basis_terms, _start_basis_at), V(_params.V)
	{
		this->delta_cdw_up = std::abs(U - V) * 0.5;
		this->delta_sc = std::abs(U + V) * 0.5;
		if (V > 0) {
			this->delta_sc *= 0.25;
		}
		else if (V < 0) {
			this->delta_cdw_up *= 0.25;
		}
		this->delta_cdw_down = (U > 0) ? -this->delta_cdw_up : this->delta_cdw_up;

		this->delta_eta = 0;// std::abs(U - V) * 0.2;
		this->delta_occupation = V * 0.1;

		this->hamilton = Matrix_L::Zero(4, 4);
	}

	Model::data_set UsingBroyden::computePhases(const bool print/*=false*/)
	{
		Matrix_L rho = Matrix_L::Zero(4, 4);
		Matrix_L trafo = Matrix_L::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<Matrix_L> solver;

		auto lambda_func = [&](const ParameterVector& x, ParameterVector& F) {
			delta_cdw_up = x(0);
			delta_cdw_down = x(1);
			delta_sc = x(2);
			delta_eta = x(3);
			delta_occupation = x(4);

			F = ParameterVector::Zero();
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * L_PI) / (Constants::K_DISCRETIZATION);
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					double_prec k_y = (l * L_PI) / (Constants::K_DISCRETIZATION);

					fillHamiltonian(k_x, k_y);
					solver.compute(hamilton);
					rho.fill(0);
					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * (solver.eigenvectors().transpose());

					F(0) -= rho(2, 3);
					F(1) += rho(0, 1);
					F(2) += rho(0, 2);
					F(3) += rho(0, 3);
					F(4) += 0.5 * cos(k_x) * (rho(0, 0) + 1 - rho(2, 2));
				}
			}
			setParameters(F);
			F -= x;
		};
		std::function<void(const ParameterVector&, ParameterVector&)> func = lambda_func;
		ParameterVector x0 = ParameterVector(delta_cdw_up, delta_cdw_down, delta_sc, delta_eta, delta_occupation);
		ParameterVector f0 = ParameterVector(delta_cdw_up, delta_cdw_down, delta_sc, delta_eta, delta_occupation);
		for (size_t i = 0; i < 200 && f0.squaredNorm() > 1e-15; i++)
		{
			func(x0, f0);
			if (std::abs(x0(0) + delta_cdw_up) < 1e-10) {
				delta_cdw_up = 0;
			}
			if (std::abs(x0(1) + delta_cdw_down) < 1e-10) {
				delta_cdw_down = 0;
			}
			if (std::abs(x0(2) + delta_sc) < 1e-10) {
				delta_sc = 0;
			}
			if (std::abs(x0(3) + delta_eta) < 1e-10) {
				delta_eta = 0;
			}
			x0(0) = delta_cdw_up;
			x0(1) = delta_cdw_down;
			x0(2) = delta_sc;
			x0(3) = delta_eta;
			x0(4) = delta_occupation;

			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8)
					<< delta_cdw_up << "\t" << delta_cdw_down << "\t" << delta_sc << "\t" << delta_eta 
					<< "\t" << delta_occupation << std::endl;
			}
		}

		Utility::NumericalSolver::Roots::Broyden<double_prec, 5> broyden_solver;
		if (!broyden_solver.compute(func, x0)) {
			std::cerr << "No convergence for [T, U, V] = [" << this->temperature << ", " << this->U << ", " << this->V << "]" << std::endl;
			delta_cdw_up = 0;
			delta_cdw_down = 0;
			delta_sc = 0;
			delta_eta = 0;
		}
		data_set ret = { this->delta_cdw_up, this->delta_cdw_down, this->delta_sc, this->delta_eta };

		if (print) {
			func(x0, f0);
			std::cout << "T=" << temperature << "   U=" << U << "   V=" << V << "\n";
			std::cout << std::scientific << std::setprecision(8) << "x0 = (";
			for (int i = 0; i < x0.size(); i++)
			{
				std::cout << " " << x0(i) << " ";
			}
			std::cout << ")\nf0 = (";
			for (int i = 0; i < f0.size(); i++)
			{
				std::cout << " " << f0(i) << " ";
			}
			std::cout << ")   -> |f0| = " << f0.norm() << std::endl;

			std::cout << "n cos = " << delta_occupation << std::endl;
		}

		return ret;
	}
}