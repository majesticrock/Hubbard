#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <iomanip>
#include "BasicHubbardModel.hpp"

namespace Hubbard {
	void BasicHubbardModel::fillHamiltonian(double_prec k_x, double_prec k_y)
	{
		hamilton.fill(0);

		hamilton(0, 1) = delta_cdw - delta_afm;
		hamilton(0, 2) = delta_sc;
		hamilton(0, 3) = I * delta_eta;
		hamilton(1, 2) = I * delta_eta;
		hamilton(1, 3) = delta_sc;
		hamilton(2, 3) = -delta_cdw - delta_afm;

		SpinorMatrix buffer = hamilton.adjoint();
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
		this->delta_cdw = 0.1;
		this->delta_afm = 0.1;
		this->delta_eta = 0;

		this->hamilton = SpinorMatrix::Zero(4, 4);
	}

	Hubbard::Model::data_set BasicHubbardModel::computePhases(const bool print/*=false*/)
	{
		SpinorMatrix rho = SpinorMatrix::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
		constexpr double_prec EPSILON = 1e-8;
		double_prec error = 100;

		auto lambda_func = [&](const ParameterVector& x, ParameterVector& F) {
			delta_cdw = x(0);
			delta_afm = x(1);
			delta_sc = x(2);
			delta_eta = x(3);

			complex_prec c_cdw = { 0, 0 }, c_afm = { 0, 0 }, c_sc = { 0, 0 }, c_eta = { 0, 0 };

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * M_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					double_prec k_y = (l * M_PI) / Constants::K_DISCRETIZATION;
					fillHamiltonian(k_x, k_y);
					solver.compute(hamilton);
					fillRho(rho, solver);

					c_cdw += (rho(0, 1) - rho(2, 3));
					c_afm += (rho(0, 1) + rho(2, 3));
					c_sc += rho(0, 2);
					c_eta += rho(0, 3);
				}
			}

			if (std::abs(c_cdw.imag()) > 1e-8) {
				std::cout << "cdw: " << c_cdw << std::endl;
			}
			if (std::abs(c_afm.imag()) > 1e-8) {
				std::cout << "afm: " << c_afm << std::endl;
			}
			if (std::abs(c_sc.imag()) > 1e-8) {
				std::cout << "sc: " << c_sc << std::endl;
			}
			if (std::abs(c_eta.real()) > 1e-8) {
				std::cout << "eta: " << c_eta << std::endl;
			}

			F(0) = c_cdw.real();
			F(1) = c_afm.real();
			F(2) = c_sc.real();
			F(3) = c_eta.imag();

			setParameters(F);
			F -= x;
		};

		constexpr int MAX_STEPS = 500;

		ParameterVector f0;
		f0 << delta_cdw, delta_afm, delta_sc, delta_eta;

		ParameterVector x0;
		x0 << delta_cdw, delta_afm, delta_sc, delta_eta;

		for (size_t i = 0; i < MAX_STEPS && error > EPSILON; i++)
		{
			lambda_func(x0, f0);
			error = f0.norm();

			x0(0) = delta_cdw;
			x0(1) = delta_afm;
			x0(2) = delta_sc;
			x0(3) = delta_eta;
			x0(4) = delta_occupation_up;

			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow(x0);
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U] = [" << this->temperature << ", " << this->U << "]\tConvergence at " << error << std::endl;
				delta_cdw = 0;
				delta_sc = 0;
				delta_eta = 0;
			}
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