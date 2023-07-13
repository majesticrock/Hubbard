#include "UsingBroyden.hpp"
#include "../../Utility/Roots_Broyden.hpp"

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

	void UsingBroyden::addToParameterSet(const SpinorMatrix& rho, ComplexParameterVector& F, const NumericalMomentum<2>& k_values)
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

	UsingBroyden::UsingBroyden(const ModelParameters& _params, size_t _MaxPreBroydenIterations/* = 300U*/)
		: Model2D(_params), MaxPreBroydenIterations(_MaxPreBroydenIterations)
	{
		init();
	}

	UsingBroyden::UsingBroyden(const ModelParameters& _params, const BaseAttributes& startingValues, size_t _MaxPreBroydenIterations/* = 300U*/)
		: Model2D(_params, startingValues), MaxPreBroydenIterations(_MaxPreBroydenIterations)
	{
		init();
	}

	ModelAttributes<double> UsingBroyden::computePhases(const PhaseDebuggingPolicy debugPolicy/*=PhaseDebuggingPolicy{}*/)
	{
		std::function<void(const ParameterVector&, ParameterVector&)> func = [&](const ParameterVector& x, ParameterVector& F) {
			iterationStep(x, F);
		};
		const size_t NUMBER_OF_PARAMETERS = model_attributes.size();
		ParameterVector f0{ ParameterVector::Zero(NUMBER_OF_PARAMETERS) };

		std::copy(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), f0.begin());
		ParameterVector x0 = f0;

		if (debugPolicy.printAll) {
			std::cout << "-1:\t" << std::fixed << std::setprecision(8);
			printAsRow<-1>(x0);
		}
		for (size_t i = 0U; i < MaxPreBroydenIterations && f0.squaredNorm() > 1e-15; ++i)
		{
			func(x0, f0);
			for (size_t j = 0U; j < NUMBER_OF_PARAMETERS; ++j)
			{
				if(std::abs(x0[j]) > 1e-10){
					if(std::abs((x0[j] + model_attributes[j]) / x0[j]) < 1e-12){
						// Sign flipping behaviour
						if (debugPolicy.convergenceWarning){
							std::cerr << "Sign flipper for [T U V] = [" << std::fixed << std::setprecision(8)
							<< this->temperature << " " << this->U << " " << this->V << "]" << std::endl;
						}

						//std::fill(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), 0.);
						//model_attributes.converged = false;
						//return ModelAttributes<double>(this->model_attributes);
					}
				}
			}
			
			
			std::copy(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), x0.begin());

			for (auto& x : x0)
			{
				if (std::abs(x) < 1e-14) {
					x = 0.;
				}
			}

			if (debugPolicy.printAll) {
				std::cout << i << ":  " << std::scientific << std::setprecision(4);
				printAsRow<-1>(x0);
			}
		}

		Utility::NumericalSolver::Roots::Broyden<double, -1> broyden_solver;
		if (!broyden_solver.compute(func, x0, 400)) {
			if (debugPolicy.convergenceWarning){
				std::cerr << "No convergence for [T U V] = [" << std::fixed << std::setprecision(8)
				<< this->temperature << " " << this->U << " " << this->V << "]" << std::endl;
			}

			std::fill(model_attributes.selfconsistency_values.begin(), model_attributes.selfconsistency_values.end(), 0.);
			model_attributes.converged = false;
		}
		else {
			model_attributes.converged = true;
		}

		if (debugPolicy.printAll) {
			func(x0, f0);
			std::cout << "T=" << temperature << "   U=" << U << "   V=" << V << "\n";
			std::cout << "x0 = (";
			for (const auto& x : x0)
			{
				std::cout << " " << x << " ";
			}
			std::cout << ")\nf0 = (";
			for (const auto& f : f0)
			{
				std::cout << " " << f << " ";
			}
			std::cout << ")\n -> |f0| = " << std::scientific << std::setprecision(8) << f0.norm() << std::endl;
		}
		return ModelAttributes<double>(this->model_attributes);
	}
}