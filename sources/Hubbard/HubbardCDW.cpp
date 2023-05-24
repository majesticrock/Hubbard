#define _USE_MATH_DEFINES

#include "HubbardCDW.hpp"

namespace Hubbard {
	void HubbardCDW::computeChemicalPotential()
	{
		Model::computeChemicalPotential();
		chemical_potential += 4 * V;
	}
	void HubbardCDW::fillHamiltonian(double_prec k_x, double_prec k_y)
	{
		hamilton.fill(0.0);

		hamilton(0, 1) = delta_cdw - delta_afm;// + I * (2 * xi_cdw_up_x * cos(k_x) + 2 * xi_cdw_up_y * cos(k_y));
		hamilton(0, 2) = delta_sc + I * (gamma_sc * gamma(k_x, k_y) + xi_sc * xi(k_x, k_y));
		hamilton(0, 3) = I * delta_eta;// +(2 * xi_eta_x * cos(k_x) + 2 * xi_eta_y * cos(k_y));

		hamilton(1, 2) = I * delta_eta;// -(2 * xi_eta_x * cos(k_x) + 2 * xi_eta_y * cos(k_y));
		hamilton(1, 3) = delta_sc - I * (gamma_sc * gamma(k_x, k_y) + xi_sc * xi(k_x, k_y));
		hamilton(2, 3) = -delta_cdw - delta_afm;// + I * (2 * xi_cdw_down_x * cos(k_x) + 2 * xi_cdw_down_y * cos(k_y));

		SpinorMatrix buffer = hamilton.adjoint();
		hamilton += buffer;
		double_prec eps = renormalizedEnergy_up(k_x, k_y);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = renormalizedEnergy_down(k_x, k_y);
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}
	HubbardCDW::HubbardCDW(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: Model(_params, _number_of_basis_terms, _start_basis_at), V(_params.V)
	{
		this->delta_cdw = std::abs(U - V) * 0.5 + 0.1;
		this->delta_sc = std::abs(U + V) * 0.5 + 0.1;
		if (V > 0) {
			this->delta_sc *= 0.25;
		}
		else if (V < 0) {
			this->delta_cdw *= 0.25;
		}
		this->delta_afm = (U > 0 && U > 4*V) ? -this->delta_cdw : this->delta_cdw;

		this->delta_eta = 0;//std::abs(U) * 0.2;
		this->delta_occupation_up = V * 0.1;
		this->delta_occupation_down = V * 0.2;
		this->delta_occupation_up_y = V * 0.3;
		this->delta_occupation_down_y = V * 0.4;
		this->gamma_sc = std::abs(V) * 0.3;
		this->xi_sc = std::abs(V) * 0.5;
		this->xi_eta_x      = std::abs(V) * 0.2;
		this->xi_eta_y      = std::abs(V) * 0.3;
		this->xi_cdw_up_x   = std::abs(V) * 0.4;
		this->xi_cdw_up_y   = std::abs(V) * 0.5;
		this->xi_cdw_down_x = std::abs(V) * 0.6;
		this->xi_cdw_down_y = std::abs(V) * 0.1;

		this->V_OVER_N = V / BASIS_SIZE;

		this->hamilton = SpinorMatrix::Zero(4, 4);
	}
	Model::data_set HubbardCDW::computePhases(const bool print)
	{
		std::cout << "V/N= " << V_OVER_N << std::endl;
		SpinorMatrix rho = SpinorMatrix::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
		constexpr double_prec EPSILON = 1e-10;
		double_prec error = 100;

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
			xi_eta_x = x(8);
			xi_eta_y = x(9);
			delta_occupation_up_y = x(10);
			delta_occupation_down_y = x(11);
			xi_cdw_up_x = x(12);
			xi_cdw_up_y = x(13);
			xi_cdw_down_x = x(14);
			xi_cdw_down_y = x(15);

			complex_prec c_sc = { 0, 0 }, c_eta = { 0, 0 };
			complex_prec c_gamma_sc = { 0, 0 }, c_xi_sc = { 0, 0 }, c_xi_eta_x = { 0, 0 }, c_xi_eta_y = { 0, 0 };
			complex_prec c_xi_cdw_up_x = { 0, 0 }, c_xi_cdw_up_y = { 0, 0 }, c_xi_cdw_down_x = { 0, 0 }, c_xi_cdw_down_y = { 0, 0 };

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * M_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < 0; l++)
				{
					double_prec k_y = (l * M_PI) / Constants::K_DISCRETIZATION;
					fillHamiltonian(k_x, k_y);
					solver.compute(hamilton);

					rho.fill(0);
					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * solver.eigenvectors().adjoint();

					F(0) += (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real();
					F(1) += (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real();
					c_sc       += (rho(2, 0) + rho(3, 1));
					c_gamma_sc += gamma(k_x, k_y) * (rho(2, 0) - rho(3, 1));
					c_xi_sc    += xi(k_x, k_y)    * (rho(2, 0) - rho(3, 1));
					c_eta      += (rho(0, 3) + rho(1, 2));
					F(6)       += cos(k_x) * (rho(0, 0) - rho(1, 1)).real();
					F(7)       -= cos(k_x) * (rho(2, 2) - rho(3, 3)).real();

					c_xi_eta_x += cos(k_x) * (rho(3, 0) - rho(2, 1));
					c_xi_eta_y += cos(k_y) * (rho(3, 0) - rho(2, 1));
					F(10)       += cos(k_y) * (rho(0, 0) - rho(1, 1)).real();
					F(11)       -= cos(k_y) * (rho(2, 2) - rho(3, 3)).real();

					c_xi_cdw_up_x   += cos(k_x) * (rho(1, 0) - rho(0, 1));
					c_xi_cdw_up_y   += cos(k_y) * (rho(1, 0) - rho(0, 1));
					c_xi_cdw_down_x -= cos(k_x) * (rho(2, 3) - rho(3, 2));
					c_xi_cdw_down_y -= cos(k_y) * (rho(2, 3) - rho(3, 2));
				}
			}

			{ // Checks for numerical accurarcy
				const double ERROR_MARGIN = EPSILON * BASIS_SIZE;
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
				if (std::abs(c_xi_eta_x.imag()) > ERROR_MARGIN) {
					std::cout << "xi eta x: " << c_xi_eta_x << std::endl;
				}
				if (std::abs(c_xi_eta_y.imag()) > ERROR_MARGIN) {
					std::cout << "xi eta y: " << c_xi_eta_y << std::endl;
				}
				if (std::abs(c_xi_cdw_up_x.real()) > ERROR_MARGIN) {
					std::cout << "xi cdw up x: " << c_xi_eta_x << std::endl;
				}
				if (std::abs(c_xi_cdw_up_y.real()) > ERROR_MARGIN) {
					std::cout << "xi cdw up y: " << c_xi_eta_y << std::endl;
				}
				if (std::abs(c_xi_cdw_down_x.real()) > ERROR_MARGIN) {
					std::cout << "xi cdw down x: " << c_xi_eta_x << std::endl;
				}
				if (std::abs(c_xi_cdw_down_y.real()) > ERROR_MARGIN) {
					std::cout << "xi cdw down y: " << c_xi_eta_y << std::endl;
				}
			}

			F(2) = c_sc.real();
			F(3) = c_gamma_sc.imag();
			F(4) = c_xi_sc.imag();
			F(5) = c_eta.imag();

			F(8) = c_xi_eta_x.real();
			F(9) = c_xi_eta_y.real();

			F(12) = c_xi_cdw_up_x.imag();
			F(13) = c_xi_cdw_up_y.imag();
			F(14) = c_xi_cdw_down_x.imag();
			F(15) = c_xi_cdw_down_y.imag();

			setParameters(F);
			F -= x;
		};

		constexpr int MAX_STEPS = 100;

		ParameterVector f0;
		f0 << delta_cdw, delta_afm, delta_sc, gamma_sc, xi_sc, delta_eta, delta_occupation_up, delta_occupation_down,
			xi_eta_x, xi_eta_y, delta_occupation_up_y, delta_occupation_down_y,
			xi_cdw_up_x, xi_cdw_up_y, xi_cdw_down_x, xi_cdw_down_y;

		ParameterVector x0;
		x0 << delta_cdw, delta_afm, delta_sc, gamma_sc, xi_sc, delta_eta, delta_occupation_up, delta_occupation_down,
			xi_eta_x, xi_eta_y, delta_occupation_up_y, delta_occupation_down_y,
			xi_cdw_up_x, xi_cdw_up_y, xi_cdw_down_x, xi_cdw_down_y;

		for (size_t i = 0; i < MAX_STEPS && error > EPSILON; i++)
		{
			lambda_func(x0, f0);
			error = f0.norm();

			x0(0) = delta_cdw;
			x0(1) = delta_afm;
			x0(2) = delta_sc;
			x0(3) = gamma_sc;
			x0(4) = xi_sc;
			x0(5) = delta_eta;
			x0(6) = delta_occupation_up;
			x0(7) = delta_occupation_down;
			x0(8) = xi_eta_x;
			x0(9) = xi_eta_y;
			x0(10) = delta_occupation_up_y;
			x0(11) = delta_occupation_down_y;
			x0(12) = xi_cdw_up_x;
			x0(13) = xi_cdw_up_y;
			x0(14) = xi_cdw_down_x;
			x0(15) = xi_cdw_down_y;

			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow(x0);
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U] = [" << this->temperature << ", " << this->U << "]\tConvergence at " << error << std::endl;
				delta_cdw = 0;
				delta_afm = 0;
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