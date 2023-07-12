#include "SquareTripletPairing.hpp"

constexpr size_t NUMBER_OF_PARAMETERS = 18;

namespace Hubbard::SquareLattice {
	void SquareTripletPairing::fillHamiltonianHelper(va_list args)
	{
		UNPACK_2D;
		hamilton.fill(0.0);
		const double GAMMA = gamma(k_x, k_y);
		const double XI = xi(k_x, k_y);

		SpinorMatrix diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 1) = DELTA_CDW - DELTA_AFM;// +((gamma_cdw - this->gamma_afm) * GAMMA + (xi_cdw - this->xi_afm) * XI);

		diagonalBlock(0, 2) = DELTA_SC + (GAMMA_SC * GAMMA + XI_SC * XI);
		diagonalBlock(0, 3) = DELTA_ETA;// +(this->gamma_eta * GAMMA + this->xi_eta * XI);

		diagonalBlock(1, 2) = DELTA_ETA;// - (this->gamma_eta * GAMMA + this->xi_eta * XI);
		diagonalBlock(1, 3) = DELTA_SC - (GAMMA_SC * GAMMA + XI_SC * XI);

		diagonalBlock(2, 3) = -DELTA_CDW - DELTA_AFM;// - ((gamma_cdw + this->gamma_afm) * GAMMA + (xi_cdw + this->xi_afm) * XI);

		SpinorMatrix buffer = diagonalBlock.adjoint();
		diagonalBlock += buffer;
		double eps = model_attributes.renormalizedEnergy_up(GAMMA) - 2. * model_attributes[12].real() * XI;
		diagonalBlock(0, 0) = eps;
		diagonalBlock(1, 1) = -eps;
		eps = model_attributes.renormalizedEnergy_down(GAMMA) - 2. * model_attributes[13].real() * XI;
		diagonalBlock(2, 2) = -eps;
		diagonalBlock(3, 3) = eps;

		hamilton.block<4, 4>(0, 0) = diagonalBlock;
		hamilton.block<4, 4>(4, 4) = -diagonalBlock.adjoint();

		const double TAU = tau(k_x, k_y);
		const double THETA = theta(k_x, k_y);
		diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 0) = tau_sc * TAU + theta_sc * THETA;
		diagonalBlock(1, 1) = -tau_sc * TAU + theta_sc * THETA;
		diagonalBlock(2, 2) = std::conj(tau_sc * TAU + theta_sc * THETA);
		diagonalBlock(3, 3) = -std::conj(tau_sc * TAU + theta_sc * THETA);

		hamilton.block<4, 4>(0, 4) = diagonalBlock;
		hamilton.block<4, 4>(4, 0) = diagonalBlock.adjoint();
	}

	SquareTripletPairing::SquareTripletPairing(const ModelParameters& _params)
		: HubbardCDW(_params)
	{
		SPINOR_SIZE = 8;
		hamilton = SpinorMatrix::Zero(8, 8);

		model_attributes.push_back(this->tau_sc);
		model_attributes.push_back(this->theta_sc);

		parameterCoefficients.push_back(V_OVER_N); // Tau_sc
		parameterCoefficients.push_back(V_OVER_N); // Theta_sc

		for (size_t i = 0; i < 16; i++)
		{
			model_attributes[i] = 0;
		}
		model_attributes[16] = (I + 0.5) * V;
		model_attributes[17] = (I + 0.5) * V;
	}

	ModelAttributes<double> SquareTripletPairing::computePhases(const bool print)
	{
		constexpr double EPSILON = 1e-12;
		double error = 100;
		constexpr size_t MAX_STEPS = 2500;

		ParameterVector f0{ ParameterVector::Zero(NUMBER_OF_PARAMETERS) };
		std::copy(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), f0.begin());

		ParameterVector x0 = f0;

		if (print) {
			std::cout << "-1:\t" << std::fixed << std::setprecision(8);
			printAsRow<-1>(x0);
		}
		for (size_t i = 0U; i < MAX_STEPS && error > EPSILON; ++i)
		{
			iterationStep(x0, f0);
			error = f0.norm();
			std::copy(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), x0.begin());

			if (print) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow<-1>(x0);
			}
			if (i == MAX_STEPS - 1) {
				std::cerr << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
					<< "]\tConvergence at " << error << std::endl;
				for (auto& value : model_attributes.selfconsistency_values) {
					value = 0;
				}
			}
		}

		if (std::abs(DELTA_SC.imag()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(GAMMA_SC.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(XI_SC.real()) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}
		if (std::abs(DELTA_ETA) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}

		return ModelAttributes<double>(this->model_attributes);
	}
}