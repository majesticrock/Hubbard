#include "TripletPairingIterative.hpp"

constexpr size_t NUMBER_OF_PARAMETERS = 18;

namespace Hubbard::SquareLattice {
	void TripletPairingIterative::fillHamiltonianHelper(va_list args)
	{
		UNPACK_2D;
		hamilton.fill(0.0);
		const double_prec GAMMA = gamma(k_x, k_y);
		const double_prec XI = xi(k_x, k_y);

		SpinorMatrix diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 1) = this->delta_cdw - this->delta_afm;// +((gamma_cdw - this->gamma_afm) * GAMMA + (xi_cdw - this->xi_afm) * XI);

		diagonalBlock(0, 2) = this->delta_sc + (this->gamma_sc * GAMMA + this->xi_sc * XI);
		diagonalBlock(0, 3) = this->delta_eta + (this->gamma_eta * GAMMA + this->xi_eta * XI);

		diagonalBlock(1, 2) = this->delta_eta - (this->gamma_eta * GAMMA + this->xi_eta * XI);
		diagonalBlock(1, 3) = this->delta_sc - (this->gamma_sc * GAMMA + this->xi_sc * XI);

		diagonalBlock(2, 3) = -this->delta_cdw - this->delta_afm;// - ((gamma_cdw + this->gamma_afm) * GAMMA + (xi_cdw + this->xi_afm) * XI);

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
		SPINOR_SIZE = 8;
		hamilton = SpinorMatrix::Zero(8, 8);

		parameterMapper.push_back(&(this->tau_sc));
		parameterMapper.push_back(&(this->theta_sc));

		parameterCoefficients.push_back(V_OVER_N); // Tau_sc
		parameterCoefficients.push_back(V_OVER_N); // Theta_sc

		for (size_t i = 0; i < 16; i++)
		{
			*(parameterMapper[i]) = 0;
		}
		*(parameterMapper[16]) = (I + 0.5) * V;
		*(parameterMapper[17]) = (I + 0.5) * V;
	}

	PhaseDataSet TripletPairingIterative::computePhases(const bool print)
	{
		constexpr double_prec EPSILON = 1e-12;
		double_prec error = 100;
		constexpr int MAX_STEPS = 2500;

		ParameterVector f0 = ParameterVector(NUMBER_OF_PARAMETERS);
		for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
		{
			f0(i) = *(parameterMapper[i]);
		}

		ParameterVector x0 = f0;

		if (print) {
			std::cout << "-1:\t" << std::fixed << std::setprecision(8);
			printAsRow(x0);
		}
		for (size_t i = 0; i < MAX_STEPS && error > EPSILON; i++)
		{
			iterationStep(x0, f0);
			error = f0.norm();
			for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
			{
				x0(i) = *(parameterMapper[i]);
			}
			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow(x0);
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
					<< "]\tConvergence at " << error << std::endl;
				delta_cdw = 0;
				delta_afm = 0;
				delta_sc = 0;
				delta_eta = 0;
			}
		}

		if (std::abs(delta_sc.imag()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(gamma_sc.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(xi_sc.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(delta_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(gamma_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(xi_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}

		PhaseDataSet ret;
		ret.delta_cdw = this->delta_cdw.real();
		ret.delta_afm = this->delta_afm.real();
		ret.delta_sc = this->delta_sc.real();
		ret.gamma_sc = this->gamma_sc.imag();
		ret.xi_sc = this->xi_sc.imag();
		ret.delta_eta = this->delta_eta.imag();

		return ret;
	}
}