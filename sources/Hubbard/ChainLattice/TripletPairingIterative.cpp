#include "TripletPairingIterative.hpp"

#define tau_sc xi_sc

namespace Hubbard::ChainLattice {
	void TripletPairingIterative::fillHamiltonianHelper(va_list args)
	{
		UNPACK_1D;
		hamilton.fill(0.0);
		const double_prec GAMMA = gamma(k_x);

		SpinorMatrix diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 1) = this->delta_cdw - this->delta_afm;

		diagonalBlock(0, 2) = this->delta_sc + this->gamma_sc * GAMMA ;
		diagonalBlock(0, 3) = this->delta_eta;

		diagonalBlock(1, 2) = this->delta_eta;
		diagonalBlock(1, 3) = this->delta_sc - this->gamma_sc * GAMMA ;

		diagonalBlock(2, 3) = -this->delta_cdw - this->delta_afm;

		SpinorMatrix buffer = diagonalBlock.adjoint();
		diagonalBlock += buffer;
		double_prec eps = renormalizedEnergy_up(GAMMA);
		diagonalBlock(0, 0) = eps;
		diagonalBlock(1, 1) = -eps;
		eps = renormalizedEnergy_down(GAMMA);
		diagonalBlock(2, 2) = -eps;
		diagonalBlock(3, 3) = eps;

		hamilton.block<4, 4>(0, 0) = diagonalBlock;
		hamilton.block<4, 4>(4, 4) = -diagonalBlock.adjoint();

		const double_prec TAU = tau(k_x);
		diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 0) = tau_sc * TAU;
		diagonalBlock(1, 1) = -tau_sc * TAU;
		diagonalBlock(2, 2) = std::conj(tau_sc * TAU);
		diagonalBlock(3, 3) = -std::conj(tau_sc * TAU);

		hamilton.block<4, 4>(0, 4) = diagonalBlock;
		hamilton.block<4, 4>(4, 0) = diagonalBlock.adjoint();
	}

	TripletPairingIterative::TripletPairingIterative(const ModelParameters& _params)
		: Model(_params)
	{
		SPINOR_SIZE = 8;
		hamilton = SpinorMatrix::Zero(8, 8);

		//*(parameterMapper[0])=0;
		//*(parameterMapper[1])=0;
		//*(parameterMapper[2]) = 0.3 * U;
		//*(parameterMapper[3]) = 0.1 * V;
		//*(parameterMapper[4]) = 0.3 * V;

		parameterCoefficients = {
			0.5 * U_OVER_N - 2. * V_OVER_N, // CDW
			0.5 * U_OVER_N, // AFM
			U_OVER_N, // SC
			2 * V_OVER_N, // Gamma SC
			2 * V_OVER_N, // Tau SC
			U_OVER_N, // Eta
			2 * V_OVER_N, // Occupation Up
			2 * V_OVER_N, // Occupation Down
		};
	}

	BaseModelRealAttributes TripletPairingIterative::computePhases(const bool print)
	{
		constexpr double_prec EPSILON = 1e-12;
		double_prec error = 100;
		constexpr int MAX_STEPS = 2000;
		const int NUMBER_OF_PARAMETERS = parameterMapper.size();
		
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
		if (std::abs(delta_eta) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}

		return BaseModelRealAttributes(*this);
	}
}