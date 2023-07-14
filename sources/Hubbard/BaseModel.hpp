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
	private:
		inline void init()
		{
			this->hamilton = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
			computeChemicalPotential();
		};

	protected:
		using ParameterVector = Eigen::Vector<DataType, Eigen::Dynamic>;
		using BaseAttributes = ModelAttributes<DataType>;

		ModelAttributes<DataType> model_attributes;
		// Stores the coefficients for the parameters (e.g. V/N) with the appropriate index
		std::vector<double> parameterCoefficients;

		SpinorMatrix hamilton;

		double temperature{};
		double U{};
		double V{};
		double U_OVER_N{ U / Constants::BASIS_SIZE };
		double V_OVER_N{ V / Constants::BASIS_SIZE };
		double chemical_potential{};

		size_t TOTAL_BASIS{};
		size_t SPINOR_SIZE{ 4 };

		inline virtual void computeChemicalPotential() {
			this->chemical_potential = 0.5 * U + 4 * V;
		};

		inline double fermi_dirac(double energy) const {
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
		inline void fillRho(SpinorMatrix& rho, const Eigen::SelfAdjointEigenSolver<SpinorMatrix>& solvedHamilton) const {
			rho.fill(0);
			for (int i = 0; i < rho.rows(); i++)
			{
				rho(i, i) = 1 - fermi_dirac(solvedHamilton.eigenvalues()(i));
			}
			rho = solvedHamilton.eigenvectors() * rho * solvedHamilton.eigenvectors().adjoint();
		};
		inline void multiplyParametersByCoefficients(ParameterVector& F) const {
			for (size_t i = 0U; i < F.size(); ++i)
			{
				F(i) *= parameterCoefficients[i];
			}
		};
		inline void setParameters(ParameterVector& F) {
			constexpr double new_weight = 0.5;
			for (size_t i = 0U; i < F.size(); ++i)
			{
				this->model_attributes[i] = new_weight * F(i) + (1 - new_weight) * this->model_attributes[i];
			}
		};

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

		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy debugPolicy = PhaseDebuggingPolicy{}) = 0;

		inline double getTotalGapValue() const {
			return this->model_attributes.getTotalGapValue();
		}
		inline virtual double entropyPerSite() = 0;
		inline virtual double internalEnergyPerSite() = 0;
		inline double freeEnergyPerSite() {
			return this->internalEnergyPerSite() - temperature * this->entropyPerSite();
		};

		inline std::string parametersAsTriplet() const {
			return ("[T U V] = [" + std::to_string(temperature) + " " + std::to_string(U) + " " + std::to_string(V) + "]");
		}
	};
}