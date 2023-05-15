#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <iomanip>
#include "BasicHubbardModel.hpp"

namespace Hubbard {
	void BasicHubbardModel::fillHamiltonian(double_prec k_x, double_prec k_y)
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

		hamilton(0, 0) = unperturbed_energy(k_x, k_y);
		hamilton(1, 1) = unperturbed_energy(k_x + M_PI, k_y + M_PI);
		hamilton(2, 2) = -unperturbed_energy(-k_x, -k_y);
		hamilton(3, 3) = -unperturbed_energy(-k_x + M_PI, -k_y + M_PI);
	}

	BasicHubbardModel::BasicHubbardModel(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: Model(_params, _number_of_basis_terms, _start_basis_at)
	{
		this->delta_sc = 0.1;
		this->delta_cdw_up = 0.1;
		this->delta_cdw_down = 0.1;
		this->delta_eta = 0;

		this->hamilton = Matrix_L::Zero(4, 4);
	}

	Hubbard::Model::data_set BasicHubbardModel::computePhases(const bool print/*=false*/)
	{
		Matrix_L rho = Matrix_L::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<Matrix_L> solver;
		constexpr double_prec EPSILON = 1e-8;
		double_prec sc = 0, cdw_up = 0, cdw_down = 0, eta = 0;
		double_prec old_parameters[4] = { 100, 100, 100, 100 };
		double_prec error = 100;
		double_prec error_cdw = 100;
		double_prec error_sc = 100;
		double_prec error_eta = 100;
		constexpr int MAX_STEPS = 10000;
		for (size_t i = 0; i < MAX_STEPS && error > EPSILON; i++)
		{
			sc = 0;
			cdw_up = 0;
			cdw_down = 0;
			eta = 0;

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					fillHamiltonian((k * M_PI) / Constants::K_DISCRETIZATION, (l * M_PI) / Constants::K_DISCRETIZATION);
					solver.compute(hamilton);
					rho.fill(0);

					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * (solver.eigenvectors().transpose());

					cdw_up -= rho(2, 3);
					cdw_down += rho(0, 1);
					sc += rho(0, 2);
					eta += rho(0, 3);
				}
			}

			old_parameters[0] = delta_cdw_up;
			old_parameters[1] = delta_cdw_down;
			old_parameters[2] = delta_sc;
			old_parameters[3] = delta_eta;
			setParameters(cdw_up, cdw_down, sc, eta);
			error_cdw = std::abs(delta_cdw_up - old_parameters[0]);
			error_sc = std::abs(delta_sc - old_parameters[1]);
			error_eta = std::abs(delta_eta - old_parameters[2]);
			error = error_cdw + error_sc + error_eta;

			delta_cdw_up = 0.5 * old_parameters[0] + 0.5 * delta_cdw_up;
			delta_cdw_up = 0.5 * old_parameters[1] + 0.5 * delta_cdw_down;
			delta_sc = 0.5 * old_parameters[2] + 0.5 * delta_sc;
			delta_eta = 0.5 * old_parameters[3] + 0.5 * delta_eta;

			if (print) {
				double_prec total = 0;
				for (size_t i = 0; i < 3; i++)
				{
					total += old_parameters[i] * old_parameters[i];
				}
				std::cout << i << ":\t" << std::fixed << std::setprecision(8)
					<< delta_cdw_up << "\t" << delta_cdw_down << "\t" << delta_sc << "\t" << delta_eta << "\t" << sqrt(total) << "\t" << error << std::endl;
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U] = [" << this->temperature << ", " << this->U << "]\tConvergence at " << error << std::endl;
				delta_cdw_up = 0;
				delta_cdw_down = 0;
				delta_sc = 0;
				delta_eta = 0;
			}
		}
		data_set ret;
		ret.delta_cdw_up = delta_cdw_up;
		ret.delta_cdw_down = delta_cdw_down;
		ret.delta_sc = delta_sc;
		ret.delta_eta = delta_eta;

		return ret;
	}
}