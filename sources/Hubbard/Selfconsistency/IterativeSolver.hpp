#pragma once
#include <iostream>
#include <iomanip>
#include "SolverBase.hpp"

namespace Hubbard::Selfconsistency {
	template <typename DataType>
	class IterativeSolver : public SolverBase<DataType> {
	protected:
		using ParameterVector = Eigen::Vector<DataType, Eigen::Dynamic>;

		inline bool hasSignFlippingBehaviour(const ParameterVector& x0) {
			for (size_t j = 0U; j < this->NUMBER_OF_PARAMETERS; ++j)
			{
				if (std::abs(x0[j]) > 1e-10) {
					if (std::abs((x0[j] + (*this->_attr)[j]) / x0[j]) < 1e-12) {
						return true;
					}
				}
			}
			return false;
		};
		bool procedureIterative(const PhaseDebuggingPolicy& debugPolicy, const size_t MAX_STEPS, const double EPSILON)
		{
			double error = 100;

			ParameterVector f0{ ParameterVector::Zero(this->NUMBER_OF_PARAMETERS) };
			std::copy(this->_attr->begin(), this->_attr->end(), f0.begin());
			ParameterVector x0{ f0 };

			if (debugPolicy.printAll) {
				std::cout << "-1:\t" << std::fixed << std::setprecision(8);
				printAsRow<-1>(x0);
			}

			size_t iterNum = 0U;
			while (iterNum < MAX_STEPS && error > EPSILON) {
				this->_model->iterationStep(x0, f0);
				if (hasSignFlippingBehaviour(x0)) {
					if (debugPolicy.convergenceWarning) {
						std::cerr << "Sign flipper for " << this->_model->parametersAsTriplet() << std::endl;
					}
					this->_attr->reset();
					return false;
				}

				error = f0.squaredNorm();
				std::copy(this->_attr->begin(), this->_attr->end(), x0.begin());

				if (debugPolicy.printAll) {
					std::cout << iterNum << ":\t" << std::fixed << std::setprecision(8);
					printAsRow<-1>(x0);
					std::cout << "Error:  " << error << "\n";
				}
				++iterNum;
			}
			if (iterNum >= MAX_STEPS) {
				return false;
			}
			return true;
		};
	public:
		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy& debugPolicy) override
		{
			constexpr double EPSILON = 1e-12;
			constexpr size_t MAX_STEPS = 1500;
			this->_attr->converged = true;
			if (!procedureIterative(debugPolicy, MAX_STEPS, EPSILON)) {
				if (debugPolicy.convergenceWarning) {
					std::cerr << "No convergence for " << this->_model->parametersAsTriplet() << std::endl;
				}
				this->_attr->reset();
			}

			return ModelAttributes<double>(*this->_attr);
		};

		IterativeSolver() = delete;
		IterativeSolver(BaseModel<DataType>* model_ptr, ModelAttributes<DataType>* attribute_ptr)
			: SolverBase<DataType>(model_ptr, attribute_ptr) {};
	};
}