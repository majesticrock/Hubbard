#define _NUM(momentum) (expecs[0](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))
#define _CDW(momentum) (expecs[1](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))
#define _SC(momentum) (expecs[2](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))
#define _ETA(momentum) (expecs[3](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))

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
	void UsingBroyden::fillHamiltonian(double k_x, double k_y)
	{
		hamilton.fill(0);

		//Eigen::Vector4d cos_sin;
		//cos_sin << cos(k_x), cos(k_y), sin(k_x), sin(k_y);
		hamilton(0, 1) = delta_cdw;
		hamilton(0, 2) = delta_sc;//+ V * old_sc.dot(cos_sin);
		hamilton(0, 3) = delta_eta;// - V * old_eta.dot(cos_sin);

		hamilton(1, 2) = delta_eta;// + V * old_eta.dot(cos_sin);
		hamilton(1, 3) = delta_sc;// -V * old_sc.dot(cos_sin);
		hamilton(2, 3) = -delta_cdw;

		Matrix_L buffer = hamilton.transpose();
		hamilton += buffer;
		//double eps = unperturbed_energy(k_x, k_y);
		double n_buffer = -2 * (cos(k_x) + cos(k_y)) * delta_occupation;
		hamilton(0, 0) = unperturbed_energy(k_x, k_y) + n_buffer;
		hamilton(1, 1) = unperturbed_energy(k_x + L_PI, k_y + L_PI) - n_buffer;
		hamilton(2, 2) = -unperturbed_energy(k_x, k_y) - n_buffer;
		hamilton(3, 3) = -unperturbed_energy(k_x + L_PI, k_y + L_PI) + n_buffer;
	}

	UsingBroyden::UsingBroyden(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: Model(_params, _number_of_basis_terms, _start_basis_at), V(_params.V)
	{
		this->delta_cdw = std::abs(U - V) * 0.5;
		this->delta_sc = std::abs(U + V) * 0.5;
		if (V > 0) {
			this->delta_sc *= 0.25;
		}
		else if (V < 0) {
			this->delta_cdw *= 0.25;
		}
		this->delta_eta = 0;
		this->delta_occupation = std::abs(V) * 0.1;

		this->hamilton = Matrix_L::Zero(4, 4);
	}

	Model::data_set UsingBroyden::computePhases(const bool print/*=false*/)
	{
		Matrix_L rho = Matrix_L::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<Matrix_L> solver;

		auto lambda_func = [&](const Eigen::VectorXd& x, Eigen::VectorXd& F) {
			delta_cdw = x(0);
			delta_sc = x(1);
			delta_eta = x(2);
			delta_occupation = x(3);

			F = Eigen::Vector4d::Zero();
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double k_x = (k * L_PI) / (Constants::K_DISCRETIZATION);
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					double k_y = (l * L_PI) / (Constants::K_DISCRETIZATION);
					fillHamiltonian(k_x, k_y);
					solver.compute(hamilton);
					rho.fill(0);
					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * (solver.eigenvectors().transpose());
					F(0) += rho(0, 1);
					F(1) += rho(0, 2);
					F(2) += rho(0, 3);
					F(3) += cos(k_x) * 0.5 * (rho(0, 0) - rho(2, 2));
					//cos_occupation += cos(k_x) * 0.5 * (rho(0, 0) - rho(2, 2));
				}
			}
			setParameters(F);
			F -= x;
		};
		std::function<void(const Eigen::VectorXd&, Eigen::VectorXd&)> func = lambda_func;
		Eigen::VectorXd x0 = Eigen::Vector4d(delta_cdw, delta_sc, delta_eta, delta_occupation);
		Eigen::VectorXd f0 = Eigen::Vector4d(delta_cdw, delta_sc, delta_eta, delta_occupation);
		for (size_t i = 0; i < 200 && f0.squaredNorm() > 1e-15; i++)
		{
			func(x0, f0);
			if (std::abs(x0(0) + delta_cdw) < 1e-10) {
				delta_cdw = 0;
			}
			if (std::abs(x0(1) + delta_sc) < 1e-10) {
				delta_sc = 0;
			}
			if (std::abs(x0(2) + delta_eta) < 1e-10) {
				delta_eta = 0;
			}
			x0(0) = delta_cdw;
			x0(1) = delta_sc;
			x0(2) = delta_eta;
			x0(3) = delta_occupation;

			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8)
					<< delta_cdw << "\t" << delta_sc << "\t" << delta_eta << "\t" << delta_occupation << std::endl;
			}
		}
		
		if (!Utility::Roots::Broyden::compute(func, x0)) {
			std::cerr << "No convergence for [T, U, V] = [" << this->temperature << ", " << this->U << ", " << this->V << "]" << std::endl;
			delta_cdw = 0;
			delta_sc = 0;
			delta_eta = 0;
		}
		data_set ret = { this->delta_cdw, this->delta_sc, this->delta_eta };

		if (print) {
			func(x0, f0);
			std::cout << "T=" << temperature << "   U=" << U << "   V=" << V << "\n";
			std::cout << "x0 = (";
			for (int i = 0; i < 3; i++)
			{
				std::cout << " " << x0(i) << " ";
			}
			std::cout << ")\nf0 = (";
			for (int i = 0; i < 3; i++)
			{
				std::cout << " " << f0(i) << " ";
			}
			std::cout << ")   -> |f0| = " << f0.norm() << std::endl;

			std::cout << "n cos = " << delta_occupation << std::endl;
		}

		return ret;
	}
}