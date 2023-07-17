#include "UsingBroyden.hpp"
#include "../Selfconsistency/BroydenSolver.hpp"

namespace Hubbard::SquareLattice {
	void UsingBroyden::init()
	{
		this->hamilton = SpinorMatrix::Zero(4, 4);

		parameterCoefficients = {
			0.5 * U_OVER_N - 4. * V_OVER_N, // CDW
			0.5 * U_OVER_N, // AFM
			U_OVER_N, // SC
			V_OVER_N, // Gamma SC
			V_OVER_N, // Xi SC
			U_OVER_N, // Eta
			V_OVER_N, // Occupation Up
			V_OVER_N, // Occupation Down
		};
	}
	void UsingBroyden::fillHamiltonian(const NumericalMomentum<2>& k_values)
	{
		hamilton.fill(0.);
		const double GAMMA = k_values.gamma();
		const double XI = xi(k_values);

		hamilton(0, 1) = DELTA_CDW - DELTA_AFM;;
		hamilton(0, 2) = DELTA_SC + (GAMMA_SC * GAMMA + I * XI_SC * XI);
		hamilton(0, 3) = I * DELTA_ETA;

		hamilton(1, 2) = I * DELTA_ETA;
		hamilton(1, 3) = DELTA_SC - (GAMMA_SC * GAMMA + I * XI_SC * XI);
		hamilton(2, 3) = -DELTA_CDW - DELTA_AFM;

		SpinorMatrix buffer{ hamilton.adjoint() };
		hamilton += buffer;
		double eps = model_attributes.renormalizedEnergy_up(GAMMA);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = model_attributes.renormalizedEnergy_down(GAMMA);
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}

	void UsingBroyden::addToParameterSet(ComplexParameterVector& F, const NumericalMomentum<2>& k_values)
	{
		const double GAMMA = k_values.gamma();
		const double XI = xi(k_values);

		F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
		F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
		F(2) -= (rho(0, 2) + rho(1, 3)); // SC
		F(3) -= GAMMA * (rho(0, 2) - rho(1, 3)); // Gamma SC
		F(4) -= XI * (rho(0, 2) - rho(1, 3)); // Xi SC
		F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
		F(6) -= GAMMA * (rho(0, 0) - rho(1, 1)).real(); // Gamma Occupation Up
		F(7) += GAMMA * (rho(2, 2) - rho(3, 3)).real(); // Gamma Occupation Down
	}

	UsingBroyden::UsingBroyden(const ModelParameters& _params, size_t MaxPreBroydenIterations)
		: Model2D(_params), _MaxPreBroydenIterations(MaxPreBroydenIterations)
	{
		init();
	}

	UsingBroyden::UsingBroyden(const ModelParameters& _params, const BaseAttributes& startingValues, size_t MaxPreBroydenIterations)
		: Model2D(_params, startingValues), _MaxPreBroydenIterations(MaxPreBroydenIterations)
	{
		init();
	}

	ModelAttributes<double> UsingBroyden::computePhases(const PhaseDebuggingPolicy debugPolicy)
	{
		Selfconsistency::BroydenSolver solver(this, &model_attributes, _MaxPreBroydenIterations);
		return solver.computePhases(debugPolicy);
	}
}