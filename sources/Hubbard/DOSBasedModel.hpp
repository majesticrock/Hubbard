#pragma once
#include "BaseModel.hpp"
#include "../../../FermionCommute/sources/Coefficient.hpp"
#include "DensityOfStates/BaseDOS.hpp"

namespace Hubbard {
	template <typename DataType>
	class DOSBasedModel : public BaseModel<DataType>
	{
	protected:
		std::shared_ptr<DensityOfStates::BaseDOS> dos;
		using ParameterVector = typename BaseModel<DataType>::ParameterVector;
		virtual inline void complexParametersToReal(const ComplexParameterVector& c, ParameterVector& r) const {
			// Does nothing, unless the derived class states otherwise
		};

		inline void iterationStep(const ParameterVector& x, ParameterVector& F) {
			F.fill(0);
			std::conditional_t<std::is_same_v<DataType, complex_prec>,
				ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			SpinorMatrix rho = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			for (size_t i = 0; i < this->model_attributes.size(); i++)
			{
				this->model_attributes[i] = x(i);
			}

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double gamma = k / Constants::K_DISCRETIZATION;
				this->fillHamiltonian(1, gamma);
				solver.compute(this->hamilton);
				this->fillRho(rho, solver);

				this->addToParameterSet(rho, complex_F, 1, gamma);
			}

			if constexpr (!std::is_same<DataType, complex_prec>::value) {
				complexParametersToReal(complex_F, F);
			}
			this->setParameters(F);
			F -= x;
		};
	public:
		DOSBasedModel(const ModelParameters& _params) : BaseModel<DataType>(_params) { };

		template<typename StartingValuesDataType>
		DOSBasedModel(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: BaseModel<DataType>(_params, startingValues) {};

		virtual inline double computeCoefficient(const SymbolicOperators::Coefficient& coeff, const double gamma) const {
			if (coeff.name == "\\epsilon_0") {
				return (-2 * gamma - this->chemical_potential);
			}
			if (coeff.name == "\\frac{U}{N}") {
				return this->U_OVER_N;
			}
			if (coeff.name == "\\tilde{V}") {
				return this->V_OVER_N * gamma;
			}
			throw(std::invalid_argument("Could not find the coefficient: " + coeff.name));
		};
	};
}