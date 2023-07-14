#include "SquareTripletPairing.hpp"
#include "../Selfconsistency/IterativeSolver.hpp"

constexpr size_t NUMBER_OF_PARAMETERS = 18;

namespace Hubbard::SquareLattice {
	void SquareTripletPairing::fillHamiltonian(const NumericalMomentum<2>& k_values)
	{
		hamilton.fill(0.0);
		const double GAMMA = k_values.gamma();
		const double XI = xi(k_values);

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

		const double TAU = k_values.tau();
		const double THETA = theta(k_values);
		diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 0) = tau_sc * TAU + theta_sc * THETA;
		diagonalBlock(1, 1) = -tau_sc * TAU + theta_sc * THETA;
		diagonalBlock(2, 2) = std::conj(tau_sc * TAU + theta_sc * THETA);
		diagonalBlock(3, 3) = -std::conj(tau_sc * TAU + theta_sc * THETA);

		hamilton.block<4, 4>(0, 4) = diagonalBlock;
		hamilton.block<4, 4>(4, 0) = diagonalBlock.adjoint();
	}

	void SquareTripletPairing::addToParameterSet(const SpinorMatrix& rho, ParameterVector& F, const NumericalMomentum<2>& k_values)
	{
		HubbardCDW::addToParameterSet(rho, F, k_values);

		F(16) -= k_values.tau() * rho(6, 2);
		F(17) -= theta(k_values[0], k_values[1]) * rho(6, 2);
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
	ModelAttributes<double> SquareTripletPairing::computePhases(const PhaseDebuggingPolicy debugPolicy)
	{
		Selfconsistency::IterativeSolver<std::complex<double>> solver(this, &model_attributes);
		return solver.computePhases(debugPolicy);
	}
}