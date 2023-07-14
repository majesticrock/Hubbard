#include "BroydenSolver.hpp"
#include "../../Utility/BroydensMethodEigen.hpp"

namespace Hubbard::Selfconsistency {
	ModelAttributes<double> BroydenSolver::computePhases(const PhaseDebuggingPolicy& debugPolicy)
	{
		procedureIterative(debugPolicy, _MaxPreBroydenIterations, 1e-15);

		std::function<void(const ParameterVector&, ParameterVector&)> func = [&](const ParameterVector& x, ParameterVector& F) {
			_model->iterationStep(x, F);
		};

		ParameterVector x0{ ParameterVector::Zero(NUMBER_OF_PARAMETERS) };
		std::copy(_attr->begin(), _attr->end(), x0.begin());
		Utility::NumericalSolver::Roots::BroydensMethodEigen<double, -1> broyden_solver;

		if (!broyden_solver.compute(func, x0, 400)) {
			if (debugPolicy.convergenceWarning) {
				std::cerr << "No convergence for " << _model->parametersAsTriplet() << std::endl;
			}
			_attr->reset();
		}
		else {
			_attr->converged = true;
		}

		if (debugPolicy.printAll) {
			ParameterVector f0{ ParameterVector::Zero(NUMBER_OF_PARAMETERS) };
			func(x0, f0);
			std::cout << _model->parametersAsTriplet() << "\n";
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

		return ModelAttributes<double>(*_attr);
	}
}