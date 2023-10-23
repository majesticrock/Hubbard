#pragma once
#include <iostream>
#include <iomanip>
#include "../GlobalDefinitions.hpp"
#include "../BaseModel.hpp"
#include "../../Utility/UnderlyingFloatingPoint.hpp"

namespace Hubbard::Selfconsistency {
	template <typename DataType>
	class IterativeSolver {
	protected:
		BaseModel<DataType>* _model{};
		ModelAttributes<DataType>* _attr{};
		const size_t NUMBER_OF_PARAMETERS;

		using ParameterVector = Eigen::Vector<DataType, Eigen::Dynamic>;

		inline bool hasSignFlippingBehaviour(const ParameterVector& x0) {
			for (size_t j = 0U; j < this->NUMBER_OF_PARAMETERS; ++j)
			{
				if (abs(x0[j]) > 1e-10) {
					if (abs((x0[j] + (*this->_attr)[j]) / x0[j]) < 1e-12) {
						return true;
					}
				}
			}
			return false;
		};
		bool procedureIterative(const PhaseDebuggingPolicy& debugPolicy, const size_t MAX_STEPS, const double EPSILON)
		{
			Utility::UnderlyingFloatingPoint_t<DataType> error{ 100 };

			ParameterVector f0{ ParameterVector::Zero(this->NUMBER_OF_PARAMETERS) };
			std::copy(this->_attr->begin(), this->_attr->end(), f0.begin());
			ParameterVector x0{ f0 };

			if (debugPolicy.printAll) {
				std::cout << "-1:\t" << std::scientific << std::setprecision(4);
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

				error = f0.norm();
				std::copy(this->_attr->begin(), this->_attr->end(), x0.begin());

				if (debugPolicy.printAll) {
					std::cout << iterNum << ":\t" << std::scientific << std::setprecision(4);
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
		virtual ModelAttributes<global_floating_type> computePhases(const PhaseDebuggingPolicy& debugPolicy)
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

			return ModelAttributes<global_floating_type>(*this->_attr);
		};

		virtual ~IterativeSolver() = default;
		IterativeSolver() = delete;
		IterativeSolver(BaseModel<DataType>* model_ptr, ModelAttributes<DataType>* attribute_ptr)
			: _model(model_ptr), _attr(attribute_ptr), NUMBER_OF_PARAMETERS(_attr->size()) {};
	};
}