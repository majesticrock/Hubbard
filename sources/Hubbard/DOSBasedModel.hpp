#pragma once
#include "BaseModel.hpp"
#include "../../../FermionCommute/sources/Coefficient.hpp"
#include "DensityOfStates/BaseDOS.hpp"

#define DELTA_CDW this->model_attributes[0]
#define DELTA_AFM this->model_attributes[1]
#define DELTA_SC this->model_attributes[2]
#define GAMMA_SC this->model_attributes[3]
// 4 is unused
#define DELTA_ETA this->model_attributes[5]
#define GAMMA_OCCUPATION_UP this->model_attributes[6]
#define GAMMA_OCCUPATION_DOWN this->model_attributes[7]


namespace Hubbard {
	template <typename DataType, typename DOS>
	class DOSBasedModel : public BaseModel<DataType>
	{
	protected:
		void fillHamiltonian(const double gamma){
			hamilton.fill(0.);

			hamilton(0, 1) = DELTA_CDW - DELTA_AFM;;
			hamilton(0, 2) = DELTA_SC + GAMMA_SC * gamma;
			hamilton(0, 3) = I * DELTA_ETA;

			hamilton(1, 2) = I * DELTA_ETA;
			hamilton(1, 3) = DELTA_SC - GAMMA_SC * gamma;
			hamilton(2, 3) = -DELTA_CDW - DELTA_AFM;

			SpinorMatrix buffer{ hamilton.adjoint() };
			hamilton += buffer;
			double eps = model_attributes.renormalizedEnergy_up(gamma);
			hamilton(0, 0) = eps;
			hamilton(1, 1) = -eps;
			eps = model_attributes.renormalizedEnergy_down(gamma);
			hamilton(2, 2) = -eps;
			hamilton(3, 3) = eps;
		};

		void addToParameterSet(const SpinorMatrix& rho, ComplexParameterVector& F, const double gamma, const size_t gamma_index){
			F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real() * DOS::values[gamma_index]; // CDW
			F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real() * DOS::values[gamma_index]; // AFM
			F(2) -= (rho(0, 2) + rho(1, 3)) * DOS::values[gamma_index]; // SC
			F(3) -= gamma * (rho(0, 2) - rho(1, 3)) * DOS::values[gamma_index]; // Gamma SC

			F(5) -= (rho(0, 3) + rho(1, 2)) * DOS::values[gamma_index]; // Eta
			F(6) -= gamma * (rho(0, 0) - rho(1, 1)).real() * DOS::values[gamma_index]; // Gamma Occupation Up
			F(7) += gamma * (rho(2, 2) - rho(3, 3)).real() * DOS::values[gamma_index]; // Gamma Occupation Down
		};

		using ParameterVector = typename BaseModel<DataType>::ParameterVector;

		inline void iterationStep(const ParameterVector& x, ParameterVector& F) {
			F.fill(0);
			std::conditional_t<std::is_same_v<DataType, complex_prec>,
				ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			SpinorMatrix rho { SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE) };
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			for (size_t i = 0U; i < this->model_attributes.size(); ++i)
			{
				this->model_attributes[i] = x(i);
			}

			for (int k = -Constants::BASIS_SIZE; k < Constants::BASIS_SIZE; ++k)
			{
				double gamma = k / Constants::BASIS_SIZE;
				this->fillHamiltonian(gamma);
				solver.compute(this->hamilton);
				this->fillRho(rho, solver);

				this->addToParameterSet(rho, complex_F, gamma, std::abs(k));
			}

			if constexpr (!std::is_same_v<DataType, complex_prec>) {
				complexParametersToReal(complex_F, F);
			}
			this->setParameters(F);
			F -= x;
		};
	public:
		DOSBasedModel(const ModelParameters& _params) : BaseModel<DataType>(_params) { this->model_attributes[4] = 0.; };

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