#define _USE_MATH_DEFINES

#include "HubbardCDW.hpp"

namespace Hubbard {
	void HubbardCDW::fillHamiltonian(double_prec k_x, double_prec k_y)
	{
		hamilton.fill(0.0);
		const double_prec GAMMA = gamma(k_x, k_y);
		const double_prec XI = xi(k_x, k_y);

		hamilton(0, 1) = this->c_delta_cdw - this->c_delta_afm;// +((gamma_cdw - this->gamma_afm) * GAMMA + (xi_cdw - this->xi_afm) * XI);

		hamilton(0, 2) = this->c_delta_sc + (this->c_gamma_sc * GAMMA + this->c_xi_sc * XI);
		hamilton(0, 3) = this->c_delta_eta + (this->c_gamma_eta * GAMMA + this->c_xi_eta * XI);

		hamilton(1, 2) = this->c_delta_eta - (this->c_gamma_eta * GAMMA + this->c_xi_eta * XI);
		hamilton(1, 3) = this->c_delta_sc - (this->c_gamma_sc * GAMMA + this->c_xi_sc * XI);

		hamilton(2, 3) = -this->c_delta_cdw - this->c_delta_afm;// - ((gamma_cdw + this->gamma_afm) * GAMMA + (xi_cdw + this->xi_afm) * XI);

		SpinorMatrix buffer = hamilton.adjoint();
		hamilton += buffer;
		double_prec eps = renormalizedEnergy_up(k_x, k_y);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = renormalizedEnergy_down(k_x, k_y);
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}
	HubbardCDW::HubbardCDW(const ModelParameters& _params)
		: Model(_params)
	{
		this->c_delta_cdw = (std::abs(U) + V) * 0.5 + 0.1;
		this->c_delta_sc =  (std::abs(U + std::abs(V)) * 0.3 + 0.05);
		if (V > 0) {
			this->c_delta_sc *= 0;
		}
		else if (V < 0) {
			this->c_delta_cdw *= 0;
		}
		if (U > 0) {
			this->c_delta_afm = std::abs(U - std::abs(V)) * 0.5 + 0.1;
		}
		else {
			this->c_delta_afm = 0;
		}
		this->c_delta_eta = 0;//U * 0.1;
		this->c_delta_occupation_up = V * 0.2;
		this->c_delta_occupation_down = V * 0.2;
		this->c_delta_occupation_up_y = V * 0.2;
		this->c_delta_occupation_down_y = V * 0.2;
		this->c_gamma_sc = I * V * 0.05;
		this->c_xi_sc = I * std::abs(V) * 0.1;
		this->c_gamma_cdw = 0;//I * V * 0.15;
		this->c_xi_cdw = 0;//I * V * 0.2;
		this->c_gamma_afm = 0;//I * V * 0.05;
		this->c_xi_afm = 0;//I * V * 0.04;
		this->c_gamma_eta = 0;//V * 0.01;
		this->c_xi_eta = 0;// V * 0.01;

		this->hamilton = SpinorMatrix::Zero(4, 4);

		param_mapper = {
			&(this->c_delta_cdw),
			&(this->c_delta_afm),
			&(this->c_delta_sc),
			&(this->c_gamma_sc),
			&(this->c_xi_sc),
			&(this->c_delta_eta),
			&(this->c_delta_occupation_up),
			&(this->c_delta_occupation_down),
			&(this->c_gamma_cdw),
			&(this->c_xi_cdw),
			&(this->c_gamma_afm),
			&(this->c_xi_afm),
			&(this->c_delta_occupation_up_y),
			&(this->c_delta_occupation_down_y),
			&(this->c_gamma_eta),
			&(this->c_xi_eta)
		};

		param_coefficients = {
			0.5 * U_OVER_N - 4. * V_OVER_N, // CDW
			0.5 * U_OVER_N, // AFM
			U_OVER_N, // SC
			V_OVER_N, // Gamma SC
			V_OVER_N, // Xi SC
			U_OVER_N, // Eta
			V_OVER_N, // Occupation Up
			V_OVER_N, // Occupation Down
			0.5 * V_OVER_N, // Gamma CDW
			0.5 * V_OVER_N, // Xi CDW
			0.5 * V_OVER_N, // Gamma AFM
			0.5 * V_OVER_N, // Xi AFM
			V_OVER_N, // Occupation Up y
			V_OVER_N, // Occupation Down y
			V_OVER_N, // Gamma SC
			V_OVER_N // Xi SC
		};
	}
	Model::data_set HubbardCDW::computePhases(const bool print)
	{
		computeChemicalPotential();

		SpinorMatrix rho = SpinorMatrix::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
		constexpr double_prec EPSILON = 1e-12;
		double_prec error = 100;

		auto lambda_func = [&](const ParameterVector& x, ParameterVector& F) {
			F.fill(0);
			for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
			{
				*(param_mapper[i]) = x(i);
			}

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * L_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < 0; l++)
				{
					double_prec k_y = (l * L_PI) / Constants::K_DISCRETIZATION;
					fillHamiltonian(k_x, k_y);
					solver.compute(hamilton);
					fillRho(rho, solver);

					F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
					F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
					F(2) -= (rho(0, 2) + rho(1, 3)); // SC
					F(3) -= gamma(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Gamma SC
					F(4) -= xi(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Xi SC
					F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
					F(6) -= cos(k_x) * (rho(0, 0) - rho(1, 1)).real(); // Occupation Up
					F(7) += cos(k_x) * (rho(2, 2) - rho(3, 3)).real(); // Occupation Down
					F(8) += I * gamma(k_x, k_y) * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Gamma CDW
					F(9) += I * xi(k_x, k_y) * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Xi CDW
					F(10) += I * gamma(k_x, k_y) * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Gamma AFM
					F(11) += I * xi(k_x, k_y) * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Xi AFM
					F(12) -= cos(k_y) * (rho(0, 0) - rho(1, 1)).real(); // Occupation Up y
					F(13) += cos(k_y) * (rho(2, 2) - rho(3, 3)).real(); // Occupation Down y
					F(14) -= gamma(k_x, k_y) * (rho(0, 3) - rho(1, 2)); // Gamma eta
					F(15) -= xi(k_x, k_y) * (rho(0, 3) - rho(1, 2)); // Xi eta
				}
			}

			for (size_t i = 0; i < F.size(); i++)
			{
				if (std::abs(F(i).real()) < 1e-12) {
					F(i) = { 0.,  F(i).imag() };
				}
				if (std::abs(F(i).imag()) < 1e-12) {
					F(i) = { F(i).real(), 0. };
				}
			}

			setParameters(F);
			F -= x;
		};

		constexpr int MAX_STEPS = 2500;

		ParameterVector f0;
		f0 << this->c_delta_cdw, this->c_delta_afm, this->c_delta_sc, this->c_gamma_sc, this->c_xi_sc, this->c_delta_eta,
			this->c_delta_occupation_up, this->c_delta_occupation_down,
			this->c_gamma_cdw, this->c_xi_cdw, this->c_gamma_afm, this->c_xi_afm,
			this->c_delta_occupation_up_y, this->c_delta_occupation_down_y,
			this->c_delta_eta, this->c_gamma_eta;

		ParameterVector x0;
		x0 << this->c_delta_cdw, this->c_delta_afm, this->c_delta_sc, this->c_gamma_sc, this->c_xi_sc, this->c_delta_eta,
			this->c_delta_occupation_up, this->c_delta_occupation_down,
			this->c_gamma_cdw, this->c_xi_cdw, this->c_gamma_afm, this->c_xi_afm,
			this->c_delta_occupation_up_y, this->c_delta_occupation_down_y,
			this->c_delta_eta, this->c_gamma_eta;

		if (print) {
			std::cout << "-1:\t" << std::fixed << std::setprecision(8);
			printAsRow(x0);
		}
		for (size_t i = 0; i < MAX_STEPS && error > EPSILON; i++)
		{
			lambda_func(x0, f0);
			error = f0.norm();
			for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
			{
				x0(i) = *(param_mapper[i]);
			}
			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow(x0);
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
					<< "]\tConvergence at " << error << std::endl;
				c_delta_cdw = 0;
				c_delta_afm = 0;
				c_delta_sc = 0;
				c_delta_eta = 0;
			}
		}

		if (print) {
			double_prec internal_energy = 0;
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * L_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < 0; l++)
				{
					double_prec k_y = (l * L_PI) / Constants::K_DISCRETIZATION;
					fillHamiltonian(k_x, k_y);
					solver.compute(hamilton);
					for (int i = 0; i < 4; i++) {
						internal_energy += (solver.eigenvalues()(i) < 0) ? solver.eigenvalues()(i) : 0;
					}
				}
			}
			std::cout << "Total energy:  " << internal_energy << std::endl;
		}

		if (std::abs(c_delta_sc.imag()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(c_gamma_sc.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(c_xi_sc.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(c_delta_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(c_gamma_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(c_xi_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}

		data_set ret;
		ret.delta_cdw = this->c_delta_cdw.real();
		ret.delta_afm = this->c_delta_afm.real();
		ret.delta_sc = this->c_delta_sc.real();
		ret.gamma_sc = this->c_gamma_sc.imag();
		ret.xi_sc = this->c_xi_sc.imag();
		ret.delta_eta = this->c_delta_eta.imag();

		return ret;
	}
}