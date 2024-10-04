#include "ChainTripletPairing.hpp"
#include <Utility/Selfconsistency/IterativeSolver.hpp>

namespace Hubbard::Models::ChainLattice {
	void ChainTripletPairing::init()
	{
		hamilton = SpinorMatrix::Zero(8, 8);
		rho = SpinorMatrix::Zero(8, 8);

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

	void ChainTripletPairing::fillHamiltonian(const NumericalMomentum<1>& k_x)
	{
		hamilton.fill(global_floating_type{});
		const global_floating_type GAMMA = k_x.gamma();

		SpinorMatrix diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 1) = DELTA_CDW - DELTA_AFM;

		diagonalBlock(0, 2) = DELTA_SC + GAMMA_SC * GAMMA;
		diagonalBlock(0, 3) = DELTA_ETA;

		diagonalBlock(1, 2) = DELTA_ETA;
		diagonalBlock(1, 3) = DELTA_SC - GAMMA_SC * GAMMA;

		diagonalBlock(2, 3) = -DELTA_CDW - DELTA_AFM;

		SpinorMatrix buffer = diagonalBlock.adjoint();
		diagonalBlock += buffer;
		global_floating_type eps = model_attributes.renormalizedEnergy_up(GAMMA);
		diagonalBlock(0, 0) = eps;
		diagonalBlock(1, 1) = -eps;
		eps = model_attributes.renormalizedEnergy_down(GAMMA);
		diagonalBlock(2, 2) = -eps;
		diagonalBlock(3, 3) = eps;

		hamilton.block<4, 4>(0, 0) = diagonalBlock;
		hamilton.block<4, 4>(4, 4) = -diagonalBlock.adjoint();

		const global_floating_type TAU = k_x.tau();
		diagonalBlock = SpinorMatrix::Zero(4, 4);
		diagonalBlock(0, 0) = TAU_SC * TAU;
		diagonalBlock(1, 1) = -TAU_SC * TAU;
		diagonalBlock(2, 2) = std::conj(TAU_SC * TAU);
		diagonalBlock(3, 3) = -std::conj(TAU_SC * TAU);

		hamilton.block<4, 4>(0, 4) = diagonalBlock;
		hamilton.block<4, 4>(4, 0) = diagonalBlock.adjoint();
	}

	void ChainTripletPairing::addToParameterSet(ParameterVector& F, const NumericalMomentum<1>& k_x)
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
		init();
	}

	ModelAttributes<global_floating_type> ChainTripletPairing::computePhases()
	{
		auto solver = Utility::Selfconsistency::make_iterative<complex_prec>(this, &model_attributes);
		return ModelAttributes<global_floating_type>(solver.compute(), Magnitude);
	}
}