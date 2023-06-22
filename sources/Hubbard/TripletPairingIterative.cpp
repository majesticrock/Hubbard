#include "TripletPairingIterative.hpp"

constexpr size_t NUMBER_OF_PARAMETERS = 18;

namespace Hubbard {
	void TripletPairingIterative::fillHamiltonian(double_prec k_x, double_prec k_y)
	{
		hamilton.fill(0.0);
		const double_prec GAMMA = gamma(k_x, k_y);
		const double_prec XI = xi(k_x, k_y);

		SpinorMatrix diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 1) = this->c_delta_cdw - this->c_delta_afm;// +((gamma_cdw - this->gamma_afm) * GAMMA + (xi_cdw - this->xi_afm) * XI);

		diagonalBlock(0, 2) = this->c_delta_sc + (this->c_gamma_sc * GAMMA + this->c_xi_sc * XI);
		diagonalBlock(0, 3) = this->c_delta_eta + (this->c_gamma_eta * GAMMA + this->c_xi_eta * XI);

		diagonalBlock(1, 2) = this->c_delta_eta - (this->c_gamma_eta * GAMMA + this->c_xi_eta * XI);
		diagonalBlock(1, 3) = this->c_delta_sc - (this->c_gamma_sc * GAMMA + this->c_xi_sc * XI);

		diagonalBlock(2, 3) = -this->c_delta_cdw - this->c_delta_afm;// - ((gamma_cdw + this->gamma_afm) * GAMMA + (xi_cdw + this->xi_afm) * XI);

		SpinorMatrix buffer = diagonalBlock.adjoint();
		diagonalBlock += buffer;
		double_prec eps = renormalizedEnergy_up(k_x, k_y);
		diagonalBlock(0, 0) = eps;
		diagonalBlock(1, 1) = -eps;
		eps = renormalizedEnergy_down(k_x, k_y);
		diagonalBlock(2, 2) = -eps;
		diagonalBlock(3, 3) = eps;

		hamilton.block<4, 4>(0, 0) = diagonalBlock;
		hamilton.block<4, 4>(4, 4) = -diagonalBlock.adjoint();

		const double_prec TAU = tau(k_x, k_y);
		const double_prec THETA = theta(k_x, k_y);
		diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 0) = tau_sc * TAU * theta_sc * THETA;
		diagonalBlock(1, 1) = -tau_sc * TAU * theta_sc * THETA;
		diagonalBlock(2, 2) = std::conj(tau_sc * TAU * theta_sc * THETA);
		diagonalBlock(3, 3) = -std::conj(tau_sc * TAU * theta_sc * THETA);

		hamilton.block<4, 4>(0, 4) = diagonalBlock;
		hamilton.block<4, 4>(4, 0) = diagonalBlock.adjoint();
	}

	TripletPairingIterative::TripletPairingIterative(const ModelParameters& _params)
		: HubbardCDW(_params)
	{
		hamilton = SpinorMatrix::Zero(8, 8);

		param_mapper.push_back(&(this->tau_sc));
		param_mapper.push_back(&(this->theta_sc));

		param_coefficients.push_back(V_OVER_N); // Tau_sc
		param_coefficients.push_back(V_OVER_N); // Theta_sc

		for (size_t i = 0; i < 16; i++)
		{
			*(param_mapper[i]) = 0;
		}
		*(param_mapper[16]) = (I + 0.5) * V;
		*(param_mapper[17]) = (I + 0.5) * V;
	}

	Model::data_set TripletPairingIterative::computePhases(const bool print)
	{
		SpinorMatrix rho = SpinorMatrix::Zero(8, 8);
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

					TripletPairingIterative::addToParameterSet(rho, F, k_x, k_y);
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
		for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
		{
			f0(i) = *(param_mapper[i]);
		}

		ParameterVector x0 = f0;

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