#include "ChainTripletPairing.hpp"

namespace Hubbard::ChainLattice {
	void ChainTripletPairing::fillHamiltonian(const NumericalMomentum<1>& k_x)
	{
		hamilton.fill(0.0);
		const double GAMMA = k_x.gamma();

		SpinorMatrix diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 1) = DELTA_CDW - DELTA_AFM;

		diagonalBlock(0, 2) = DELTA_SC + GAMMA_SC * GAMMA;
		diagonalBlock(0, 3) = DELTA_ETA;

		diagonalBlock(1, 2) = DELTA_ETA;
		diagonalBlock(1, 3) = DELTA_SC - GAMMA_SC * GAMMA;

		diagonalBlock(2, 3) = -DELTA_CDW - DELTA_AFM;

		SpinorMatrix buffer = diagonalBlock.adjoint();
		diagonalBlock += buffer;
		double eps = model_attributes.renormalizedEnergy_up(GAMMA);
		diagonalBlock(0, 0) = eps;
		diagonalBlock(1, 1) = -eps;
		eps = model_attributes.renormalizedEnergy_down(GAMMA);
		diagonalBlock(2, 2) = -eps;
		diagonalBlock(3, 3) = eps;

		hamilton.block<4, 4>(0, 0) = diagonalBlock;
		hamilton.block<4, 4>(4, 4) = -diagonalBlock.adjoint();

		const double TAU = k_x.tau();
		diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 0) = TAU_SC * TAU;
		diagonalBlock(1, 1) = -TAU_SC * TAU;
		diagonalBlock(2, 2) = std::conj(TAU_SC * TAU);
		diagonalBlock(3, 3) = -std::conj(TAU_SC * TAU);

		hamilton.block<4, 4>(0, 4) = diagonalBlock;
		hamilton.block<4, 4>(4, 0) = diagonalBlock.adjoint();
	}

	void ChainTripletPairing::addToParameterSet(const SpinorMatrix& rho, ParameterVector& F, const NumericalMomentum<1>& k_x)
	{
		F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
		F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
		F(2) -= (rho(0, 2) + rho(1, 3)); // SC
		F(3) -= k_x.gamma() * (rho(0, 2) - rho(1, 3)); // Gamma SC
		F(4) -= k_x.tau() * rho(6, 2); // Tau SC
		F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
		F(6) -= k_x.gamma() * (rho(0, 0) - rho(1, 1)).real(); // Gamma Occupation Up
		F(7) += k_x.gamma() * (rho(2, 2) - rho(3, 3)).real(); // Gamma Occupation Down
	}

	ChainTripletPairing::ChainTripletPairing(const ModelParameters& _params)
		: Model1D(_params)
	{
		SPINOR_SIZE = 8;
		hamilton = SpinorMatrix::Zero(8, 8);

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

	ModelAttributes<double> ChainTripletPairing::computePhases(const PhaseDebuggingPolicy debugPolicy/*=PhaseDebuggingPolicy{}*/)
	{
		constexpr double EPSILON = 1e-12;
		double error = 100;
		constexpr size_t MAX_STEPS = 2000;
		const size_t NUMBER_OF_PARAMETERS = model_attributes.size();

		ParameterVector f0{ ParameterVector::Zero(NUMBER_OF_PARAMETERS) };
		std::copy(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), f0.begin());

		ParameterVector x0 = f0;

		if (debugPolicy.printAll) {
			std::cout << "-1:\t" << std::fixed << std::setprecision(8);
			printAsRow<-1>(x0);
		}
		for (size_t i = 0U; i < MAX_STEPS && error > EPSILON; ++i)
		{
			iterationStep(x0, f0);
			error = f0.norm();
			std::copy(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), x0.begin());

			if (debugPolicy.printAll) {
				std::cout << i << ":\t" << std::fixed << std::setprecision(8);
				printAsRow<-1>(x0);
			}
			if (i == MAX_STEPS - 1) {
				if (debugPolicy.convergenceWarning){
					std::cerr << "No convergence for [T U V] = [" << std::fixed << std::setprecision(8)
					<< this->temperature << " " << this->U << " " << this->V << "]" << std::endl;
				}

				std::fill(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), 0.);
				model_attributes.converged = false;
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
		if (std::abs(DELTA_ETA) > 1e-8) {
			std::cout << "[T, U, V] = [" << this->temperature << ", " << this->U << "," << this->V
				<< "]" << std::endl;
		}

		return ModelAttributes<double>(this->model_attributes);
	}
}