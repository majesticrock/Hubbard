#include "Model.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <fstream>
#include <algorithm>
#include "../Utility/Lanczos.hpp"

// Both methods yield precisely the same data!
#define _PSEUDO_INVERSE

namespace Hubbard {
	constexpr double_prec SQRT_SALT = 1e-5;
	constexpr double_prec SALT = SQRT_SALT * SQRT_SALT;
	constexpr double_prec ERROR_MARGIN = 1e-10;

	void Model::initializeParameters()
	{
		this->U_OVER_N = U / Constants::BASIS_SIZE;
		this->V_OVER_N = V / Constants::BASIS_SIZE;
		this->delta_cdw = 0.1;
		this->delta_afm = (U > 0) ? -this->delta_cdw : this->delta_cdw;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
		computeChemicalPotential();
	}

	Model::Model(const ModelParameters& _params)
		: temperature(_params.temperature), U(_params.U), V(_params.V)
	{
		initializeParameters();
	}

	void Model::getEnergies(std::vector<std::vector<double>>& reciever, double_prec direction)
	{
		reciever.reserve(2 * Constants::K_DISCRETIZATION);
		Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
		double_prec k_val = 0;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			k_val = k * L_PI / Constants::K_DISCRETIZATION;
			fillHamiltonian(cos(L_PI * direction) * k_val, sin(L_PI * direction) * k_val);
			solver.compute(hamilton, false);
			reciever.push_back(std::vector<double>(solver.eigenvalues().data(), solver.eigenvalues().data() + solver.eigenvalues().size()));
		}
	}

	void Model::getAllEnergies(std::vector<std::vector<double>>& reciever)
	{
		reciever.resize(4 * Constants::K_DISCRETIZATION, std::vector<double>(2 * Constants::K_DISCRETIZATION));
		Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
		double_prec k_val = 0;
		double_prec l_val = 0;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			k_val = k * L_PI / Constants::K_DISCRETIZATION;
			for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
			{
				l_val = l * L_PI / Constants::K_DISCRETIZATION;
				fillHamiltonian(k_val, l_val);
				solver.compute(hamilton, false);
				reciever[k + Constants::K_DISCRETIZATION][l + Constants::K_DISCRETIZATION] = solver.eigenvalues()(0);

				for (int i = 1; i < 4; i++)
				{
					if (std::abs(solver.eigenvalues()(0) - solver.eigenvalues()(i)) > 1e-8) {
						reciever[k + 3 * Constants::K_DISCRETIZATION][l + Constants::K_DISCRETIZATION] = solver.eigenvalues()(i);
						break;
					}
				}
			}
		}
	}

	void Model::computeExpectationValues(std::vector<Matrix_L>& expecs, std::vector<double>& sum_of_all)
	{
		SpinorMatrix rho = SpinorMatrix::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

		expecs = std::vector<Matrix_L>(5, Matrix_L::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
		sum_of_all = std::vector<double>(5, 0.0);

		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
			{
				fillHamiltonian((k * L_PI) / Constants::K_DISCRETIZATION, (l * L_PI) / Constants::K_DISCRETIZATION);
				solver.compute(hamilton);
				rho.fill(0);
				for (int i = 0; i < 4; i++)
				{
					rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
				}
				rho = solver.eigenvectors() * rho * (solver.eigenvectors().adjoint());

				expecs[0](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = 1 - rho(0, 0).real();
				expecs[1](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = -rho(1, 0).real();
				expecs[2](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = -rho(0, 2).real();
				expecs[3](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = -rho(0, 3).real();
				expecs[4](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = rho(2, 3).real();
				for (int idx = 0; idx < 4; idx++)
				{
					sum_of_all[idx] += expecs[idx](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION);

					if (std::abs(rho(idx, 0).imag()) > SALT) {
						std::cerr << "Expectation values are complex!  " << rho(idx, 0) << std::endl;
						throw;
					}
				}

				if (std::abs(rho(3, 0)) > SALT) {
					std::cerr << "Warning: <eta> does not vanish! " << rho(3, 0) << std::endl;
				}
			}
		}
	}
}