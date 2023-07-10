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
	void UsingBroyden::fillHamiltonianHelper(va_list args)
	{
		UNPACK_2D;
		hamilton.fill(0);
		const double GAMMA = gamma(k_x, k_y);
		const double XI = xi(k_x, k_y);

		hamilton(0, 1) = DELTA_CDW - DELTA_AFM;;
		hamilton(0, 2) = DELTA_SC + (GAMMA_SC * GAMMA + I * XI_SC * XI);
		hamilton(0, 3) = I * DELTA_ETA;

		hamilton(1, 2) = I * DELTA_ETA;
		hamilton(1, 3) = DELTA_SC - (GAMMA_SC * GAMMA + I * XI_SC * XI);
		hamilton(2, 3) = -DELTA_CDW - DELTA_AFM;

		SpinorMatrix buffer = hamilton.adjoint();
		hamilton += buffer;
		double eps = model_attributes.renormalizedEnergy_up(GAMMA);
		hamilton(0, 0) = eps;
		hamilton(1, 1) = -eps;
		eps = model_attributes.renormalizedEnergy_down(GAMMA);
		hamilton(2, 2) = -eps;
		hamilton(3, 3) = eps;
	}

	UsingBroyden::UsingBroyden(const ModelParameters& _params, int _MaxPreBroydenIterations/* = 300*/)
		: Model2D(_params), MaxPreBroydenIterations(_MaxPreBroydenIterations)
	{
		init();
	}

	UsingBroyden::UsingBroyden(const ModelParameters& _params, const BaseAttributes& startingValues, int _MaxPreBroydenIterations/* = 300*/)
		: Model2D(_params, startingValues), MaxPreBroydenIterations(_MaxPreBroydenIterations)
	{
		init();
	}

	ModelAttributes<double> UsingBroyden::computePhases(const bool print/*=false*/)
	{
		std::function<void(const ParameterVector&, ParameterVector&)> func = [&](const ParameterVector& x, ParameterVector& F) {
			iterationStep(x, F);
		};
		const size_t NUMBER_OF_PARAMETERS = model_attributes.size();
		ParameterVector f0 = ParameterVector::Zero(NUMBER_OF_PARAMETERS);
		for (size_t i = 0U; i < NUMBER_OF_PARAMETERS; ++i)
		{
			f0(i) = model_attributes[i];
		}
		ParameterVector x0 = f0;

		for (size_t i = 0U; i < MaxPreBroydenIterations && f0.squaredNorm() > 1e-15; ++i)
		{
			func(x0, f0);
			for (size_t j = 0U; j < NUMBER_OF_PARAMETERS; ++j)
			{
				x0(j) = model_attributes[j];
			}

			for (size_t i = 0; i < x0.size(); i++)
			{
				if (std::abs(x0(i)) < 1e-14) {
					x0(i) = 0;
				}
			}

			if (print) {
				std::cout << i << ":  " << std::scientific << std::setprecision(4);
				printAsRow<-1>(x0);
			}
		}

		Utility::NumericalSolver::Roots::Broyden<double, -1> broyden_solver;
		if (!broyden_solver.compute(func, x0, 400)) {
			std::cerr << "No convergence for [T U V] = [" << std::fixed << std::setprecision(8)
				<< this->temperature << " " << this->U << " " << this->V << "]" << std::endl;
			for (auto& value : model_attributes.selfconsistency_values) {
				value = 0.;
			}
			model_attributes.converged = false;
		}
		else {
			model_attributes.converged = true;
		}

		if (print) {
			func(x0, f0);
			std::cout << "T=" << temperature << "   U=" << U << "   V=" << V << "\n";
			std::cout << "x0 = (";
			for (int i = 0; i < x0.size(); i++)
			{
				std::cout << " " << x0(i) << " ";
			}
			std::cout << ")\nf0 = (";
			for (int i = 0; i < f0.size(); i++)
			{
				std::cout << " " << f0(i) << " ";
			}
			std::cout << ")\n -> |f0| = " << std::scientific << std::setprecision(8) << f0.norm() << std::endl;
		}
		return ModelAttributes<double>(this->model_attributes);
	}
}