#pragma once
#include "GlobalDefinitions.hpp"
#include "Constants.hpp"
#include "ModelAttributes.hpp"
#include <string>
#include <memory>

namespace Hubbard {
	template <typename DataType>
	class BaseModel
	{
	public:
		using ParameterVector = Eigen::Vector<DataType, Eigen::Dynamic>;
		using BaseAttributes = ModelAttributes<DataType>;
		using HamiltonSolver = Eigen::SelfAdjointEigenSolver<SpinorMatrix>;
	private:
		inline void init()
		{
			this->hamilton = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
			this->rho = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);

			computeChemicalPotential();
			this->parameterCoefficients = {
				0.5 * this->U_OVER_N - 4. * this->V_OVER_N, // CDW
				0.5 * this->U_OVER_N, // AFM
				this->U_OVER_N, // SC
				this->V_OVER_N, // Gamma SC
				this->V_OVER_N, // Xi SC
				this->U_OVER_N, // Eta
				this->V_OVER_N, // Occupation Up
				this->V_OVER_N, // Occupation Down
			};
		};
		inline void multiplyParametersByCoefficients(ParameterVector& F) const {
			for (size_t i = 0U; i < F.size(); ++i)
			{
				F(i) *= parameterCoefficients[i];
			}
		};
		inline void setParameters(ParameterVector& F) {
			constexpr global_floating_type new_weight{ 0.5 };
			for (size_t i = 0U; i < F.size(); ++i)
			{
				this->model_attributes[i] = new_weight * F(i) + (1 - new_weight) * this->model_attributes[i];
			}
		};

	protected:
		ModelAttributes<DataType> model_attributes;
		// Stores the coefficients for the parameters (e.g. V/N) with the appropriate index
		std::vector<global_floating_type> parameterCoefficients;

		SpinorMatrix hamilton;
		SpinorMatrix rho;
		HamiltonSolver hamilton_solver;

		global_floating_type temperature{};
		global_floating_type U{};
		global_floating_type V{};
		global_floating_type U_OVER_N{ U / Constants::BASIS_SIZE };
		global_floating_type V_OVER_N{ V / Constants::BASIS_SIZE };
		global_floating_type chemical_potential{};

		size_t TOTAL_BASIS{};
		size_t SPINOR_SIZE{ 4 };

		inline virtual void computeChemicalPotential() {
			this->chemical_potential = 0.5 * U + 4 * V;
		};

		inline global_floating_type fermi_dirac(global_floating_type energy) const {
			if (temperature > 1e-8) {
				return (1. / (1. + exp(energy / temperature)));
			}
			else {
				if (std::abs(energy) < 1e-12) {
					return 0.5;
				}
				return ((energy > 0) ? 0 : 1);
			}
		};
		inline void fillRho() {
			this->hamilton_solver.compute(this->hamilton);
			rho.fill(0);
			for (int i = 0; i < rho.rows(); i++)
			{
				rho(i, i) = 1 - fermi_dirac(hamilton_solver.eigenvalues()(i));
			}
			rho = hamilton_solver.eigenvectors() * rho * hamilton_solver.eigenvectors().adjoint();
		};
		inline void applyIteration(ParameterVector& F) {
			this->multiplyParametersByCoefficients(F);
			// Numerical noise correction
			for (auto& value : F) {
				if (std::abs(value) < 1e-12) value = 0.;
			}
			this->setParameters(F);
		}
	public:
		explicit BaseModel(const ModelParameters& params)
			: model_attributes(params), temperature(params.temperature), U(params.U), V(params.V)
		{
			init();
		};

		template<typename StartingValuesDataType>
		BaseModel(const ModelParameters& _params, const ModelAttributes< StartingValuesDataType >& startingValues)
			: model_attributes(startingValues), temperature(_params.temperature), U(_params.U), V(_params.V)
		{
			init();
		};
		virtual ~BaseModel() = default;

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) = 0;

		virtual ModelAttributes<global_floating_type> computePhases(const PhaseDebuggingPolicy debugPolicy = PhaseDebuggingPolicy{}) = 0;

		inline auto getTotalGapValue() const {
			return this->model_attributes.getTotalGapValue();
		}
		inline virtual global_floating_type entropyPerSite() = 0;
		inline virtual global_floating_type internalEnergyPerSite() = 0;
		inline global_floating_type freeEnergyPerSite() {
			return this->internalEnergyPerSite() - temperature * this->entropyPerSite();
		};

		inline std::string parametersAsTriplet() const {
			return ("[T U V] = [" + std::to_string(temperature) + " " + std::to_string(U) + " " + std::to_string(V) + "]");
		}
	};
}