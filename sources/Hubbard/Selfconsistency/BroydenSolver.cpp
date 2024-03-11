#include "BroydenSolver.hpp"
#include "../../../../Utility/sources/BroydensMethodEigen.hpp"

namespace Hubbard::Selfconsistency {
	const ModelAttributes<global_floating_type>& BroydenSolver::computePhases(const PhaseDebuggingPolicy& debugPolicy)
	{
		procedureIterative(debugPolicy, _MaxPreBroydenIterations, DEFAULT_PRECISION);

		std::function<void(const ParameterVector&, ParameterVector&)> func = [&](const ParameterVector& x, ParameterVector& F) {
			_model->iterationStep(x, F);
			};

		ParameterVector x0{ ParameterVector::Zero(NUMBER_OF_PARAMETERS) };
		std::copy(_attr->begin(), _attr->end(), x0.begin());
		Utility::NumericalSolver::Roots::BroydensMethodEigen<global_floating_type, -1> broyden_solver;

		if (!broyden_solver.compute(func, x0, 400)) {
			if (debugPolicy.convergenceWarning) {
				std::cerr << std::fixed << std::setprecision(8) << "No convergence for " << _model->parametersAsTriplet() << std::endl;
			}
			_attr->reset();
		}
		else {
			_attr->converged = true;
		}

		if (debugPolicy.printAll) {
			ParameterVector f0{ ParameterVector::Zero(NUMBER_OF_PARAMETERS) };
			_model->iterationStep(x0, f0);
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

		return *this->_attr;
	}
}