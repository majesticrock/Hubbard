#define _USE_MATH_DEFINES

#include "HubbardCDW.hpp"

namespace Hubbard {
	void HubbardCDW::computeChemicalPotential()
	{
		Model::computeChemicalPotential();
		chemical_potential += 0.5 * V;
	}
	void HubbardCDW::fillHamiltonian(double_prec k_x, double_prec k_y)
	{
		complex_h.fill(0);

		complex_h(0, 1) = delta_cdw_up;
		complex_h(0, 2) = delta_sc - I * (2 * xi_sc_x * cos(k_x) + 2 * xi_sc_y * cos(k_y));
		complex_h(0, 3) = I * delta_eta - 2 * xi_eta_x * cos(k_x) - 2 * xi_eta_y * cos(k_y);

		complex_h(1, 2) = I * delta_eta + 2 * xi_eta_x * cos(k_x) - 2 * xi_eta_y * cos(k_y);
		complex_h(1, 3) = delta_sc + I * (2 * xi_sc_x * cos(k_x) + 2 * xi_sc_y * cos(k_y));
		complex_h(2, 3) = -delta_cdw_down;

		Matrix_4cL buffer = complex_h.adjoint();
		complex_h += buffer;
		double_prec eps = renormalizedEnergy_up(k_x, k_y);
		complex_h(0, 0) = eps;
		complex_h(1, 1) = -eps;
		eps = renormalizedEnergy_down(k_x, k_y);
		complex_h(2, 2) = -eps;
		complex_h(3, 3) = eps;
	}
	HubbardCDW::HubbardCDW(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: Model(_params, _number_of_basis_terms, _start_basis_at), V(_params.V)
	{
		this->delta_cdw_up = std::abs(U - V) * 0.5 + 0.1;
		this->delta_sc = std::abs(U + V) * 0.5 + 0.1;
		if (V > 0) {
			this->delta_sc *= 0.25;
		}
		else if (V < 0) {
			this->delta_cdw_up *= 0.25;
		}
		this->delta_cdw_down = ((U - V) >= 0) ? -this->delta_cdw_up : this->delta_cdw_up;

		this->delta_eta = std::abs(U - V) * 0.2;
		this->delta_occupation_up = V * 0.1;
		this->delta_occupation_down = V * 0.1;
		this->xi_sc_x = this->delta_sc * 0.1;
		this->xi_sc_y = -this->delta_sc * 0.2;
		this->xi_eta_x = this->delta_sc * 0.1;
		this->xi_eta_y = -this->delta_sc * 0.2;

		this->hamilton = Matrix_L::Zero(4, 4);
	}
	Model::data_set HubbardCDW::computePhases(const bool print)
	{
		Matrix_4cL rho = Matrix_4cL::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<Matrix_4cL> solver;
		constexpr double_prec EPSILON = 1e-8;
		double_prec error = 100;

		auto lambda_func = [&](const ParameterVector& x, ParameterVector& F) {
			delta_cdw_up = x(0);
			delta_cdw_down = x(1);
			delta_sc = x(2);
			delta_eta = x(3);
			delta_occupation_up = x(4);
			delta_occupation_down = x(5);
			xi_sc_x = x(6);
			xi_sc_y = x(7);
			xi_eta_x = x(8);
			xi_eta_y = x(9);

			complex_prec c_cdw_up = { 0, 0 }, c_cdw_down = { 0, 0 }, c_sc = { 0, 0 }, c_eta = { 0, 0 };
			complex_prec c_xi_sc_x = { 0,0 }, c_xi_sc_y = { 0,0 }, c_xi_eta_x = { 0,0 }, c_xi_eta_y = { 0,0 };

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * M_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					double_prec k_y = (l * M_PI) / Constants::K_DISCRETIZATION;
					fillHamiltonian(k_x, k_y);
					solver.compute(complex_h);

					rho.fill(0);
					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * solver.eigenvectors().adjoint();

					c_cdw_up -= rho(2, 3);
					c_cdw_down += rho(0, 1);
					c_sc += rho(0, 2);
					c_eta += rho(0, 3);
					F(4) += cos(k_x) * rho(0, 0).real();
					F(5) += cos(k_x) * (1 - rho(2, 2).real());
					c_xi_sc_x += cos(k_x) * rho(0, 2);
					c_xi_eta_x += cos(k_x) * rho(0, 3);
					c_xi_sc_y += cos(k_y) * rho(0, 2);
					c_xi_eta_y += cos(k_y) * rho(0, 3);
				}
			}

			if (std::abs(c_cdw_up.imag()) > 1e-8) {
				std::cout << "cdw_up: " << c_cdw_up << std::endl;
			}
			if (std::abs(c_cdw_down.imag()) > 1e-8) {
				std::cout << "cdw_down: " << c_cdw_down << std::endl;
			}
			if (std::abs(c_sc.imag()) > 1e-8) {
				std::cout << "sc: " << c_sc << std::endl;
			}
			if (std::abs(c_eta.real()) > 1e-8) {
				std::cout << "eta: " << c_eta << std::endl;
			}
			if (std::abs(c_xi_sc_x.real()) > 1e-8) {
				std::cout << "xi sc x: " << c_xi_sc_x << std::endl;
			}
			if (std::abs(c_xi_sc_y.real()) > 1e-8) {
				std::cout << "xi sc y: " << c_xi_sc_y << std::endl;
			}
			if (std::abs(c_xi_eta_x.imag()) > 1e-8) {
				std::cout << "xi eta x: " << c_xi_eta_x << std::endl;
			}
			if (std::abs(c_xi_eta_y.imag()) > 1e-8) {
				std::cout << "xi eta y: " << c_xi_eta_y << std::endl;
			}

			F(0) = c_cdw_up.real();
			F(1) = c_cdw_down.real();
			F(2) = c_sc.real();
			F(3) = c_eta.imag();
			F(6) = c_xi_sc_x.imag();
			F(7) = c_xi_sc_y.imag();
			F(8) = c_xi_eta_x.real();
			F(9) = c_xi_eta_y.real();
			setParameters(F);
			F -= x;
		};

		constexpr int MAX_STEPS = 1000;

		ParameterVector f0;
		f0 << delta_cdw_up, delta_cdw_down, delta_sc, delta_eta, delta_occupation_up, delta_occupation_down, 
			xi_sc_x, xi_sc_y, xi_eta_x, xi_eta_y;
		ParameterVector x0;
		x0 << delta_cdw_up, delta_cdw_down, delta_sc, delta_eta, delta_occupation_up, delta_occupation_down, 
			xi_sc_x, xi_sc_y, xi_eta_x, xi_eta_y;

		for (size_t i = 0; i < MAX_STEPS && error > EPSILON; i++)
		{
			lambda_func(x0, f0);
			error = f0.norm();

			x0(0) = delta_cdw_up;
			x0(1) = delta_cdw_down;
			x0(2) = delta_sc;
			x0(3) = delta_eta;
			x0(4) = delta_occupation_up;
			x0(5) = delta_occupation_down;
			x0(6) = xi_sc_x;
			x0(7) = xi_sc_y;
			x0(8) = xi_eta_x;
			x0(9) = xi_eta_y;

			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8)
					<< delta_cdw_up << "\t" << delta_cdw_down << "\t" << delta_sc << "\t" << delta_eta
					<< "\t" << delta_occupation_up << "\t" << delta_occupation_down
					<< "\t" << xi_sc_x << "\t" << xi_sc_y << "  " << xi_eta_x << "  " << xi_eta_y << std::endl;
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U] = [" << this->temperature << ", " << this->U << "]\tConvergence at " << error << std::endl;
				delta_cdw_up = 0;
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