#pragma once
#include "../BaseModel.hpp"
#include "../../../../FermionCommute/sources/Coefficient.hpp"
#include "../DensityOfStates/Square.hpp"
#include <algorithm>

#define DELTA_CDW this->model_attributes[0]
#define DELTA_AFM this->model_attributes[1]
#define DELTA_SC this->model_attributes[2]
#define GAMMA_SC this->model_attributes[3]
// 4 is unused
#define DELTA_ETA this->model_attributes[5]
#define GAMMA_OCCUPATION_UP this->model_attributes[6]
#define GAMMA_OCCUPATION_DOWN this->model_attributes[7]

namespace Hubbard {
	template <typename DataType, class DOS>
	class DOSBasedModel : public BaseModel<DataType>
	{
	private:
		void init() {
			this->model_attributes[4] = 0.;
			if (!DOS::computed) {
				DOS dos;
				dos.computeValues();
			}

			this->U_OVER_N = 2 * this->U / (2 * Constants::BASIS_SIZE - 1);
			this->V_OVER_N = 2 * this->V / (2 * Constants::BASIS_SIZE - 1);

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
	protected:
		using ParameterVector = typename BaseModel<DataType>::ParameterVector;

		void fillHamiltonian(const double gamma) {
			this->hamilton.fill(0.);

			this->hamilton(0, 1) = DELTA_CDW - DELTA_AFM;;
			this->hamilton(0, 2) = DELTA_SC + GAMMA_SC * gamma;
			this->hamilton(0, 3) = I * DELTA_ETA;

			this->hamilton(1, 2) = I * DELTA_ETA;
			this->hamilton(1, 3) = DELTA_SC - GAMMA_SC * gamma;
			this->hamilton(2, 3) = -DELTA_CDW - DELTA_AFM;

			SpinorMatrix buffer{ this->hamilton.adjoint() };
			this->hamilton += buffer;
			double eps = this->model_attributes.renormalizedEnergy_up(gamma);
			this->hamilton(0, 0) = eps;
			this->hamilton(1, 1) = -eps;
			eps = this->model_attributes.renormalizedEnergy_down(gamma);
			this->hamilton(2, 2) = -eps;
			this->hamilton(3, 3) = eps;
		};

		void addToParameterSet(const SpinorMatrix& rho, ComplexParameterVector& F, const double gamma, const double dos_value) {
			F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real() * dos_value; // CDW
			F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real() * dos_value; // AFM
			F(2) -= (rho(0, 2) + rho(1, 3)) * dos_value; // SC
			F(3) -= gamma * (rho(0, 2) - rho(1, 3)) * dos_value; // Gamma SC

			F(5) -= (rho(0, 3) + rho(1, 2)) * dos_value; // Eta
			F(6) -= gamma * (rho(0, 0) - rho(1, 1)).real() * dos_value; // Gamma Occupation Up
			F(7) += gamma * (rho(2, 2) - rho(3, 3)).real() * dos_value; // Gamma Occupation Down
		};

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override {
			F.fill(0);
			std::conditional_t<std::is_same_v<DataType, complex_prec>,
				ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			SpinorMatrix rho{ SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE) };
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			std::copy(x.begin(), x.end(), this->model_attributes.selfconsistency_values.begin());

			for (int k = -Constants::BASIS_SIZE + 1; k < Constants::BASIS_SIZE; ++k)
			{
				double gamma = (0.5 + k) * DOS::step;
				this->fillHamiltonian(gamma);
				solver.compute(this->hamilton);
				this->fillRho(rho, solver);

				this->addToParameterSet(rho, complex_F, gamma, DOS::values[std::abs(k)]);
			}

			if constexpr (!std::is_same_v<DataType, complex_prec>) {
				complexParametersToReal(complex_F, F);
			}
			this->applyIteration(F);
			F -= x;
		};
	public:
		DOSBasedModel(const ModelParameters& _params) : BaseModel<DataType>(_params)
		{
			init();
		};

		template<typename StartingValuesDataType>
		DOSBasedModel(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: BaseModel<DataType>(_params, startingValues)
		{
			init();
		};

		inline virtual double entropyPerSite() override {
			double entropy = 0;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			for (int k = -Constants::BASIS_SIZE + 1; k < Constants::BASIS_SIZE; k++)
			{
				double gamma = (0.5 + k) * DOS::step;
				this->fillHamiltonian(gamma);
				solver.compute(this->hamilton, false);

				entropy += std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), double{},
					[this, k](double current, double toAdd) {
						auto occ = BaseModel<DataType>::fermi_dirac(toAdd);
						// Let's just not take the ln of 0. Negative numbers cannot be reached (because math...)
						return (occ > 1e-12 ? current - occ * std::log(occ) : current) * DOS::values[std::abs(k)];
					});
			}
			return 2 * entropy / (2 * Constants::BASIS_SIZE - 1);
		};

		inline virtual double internalEnergyPerSite() override {
			double energy = 0;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			for (int k = -Constants::BASIS_SIZE + 1; k < Constants::BASIS_SIZE; k++)
			{
				double gamma = (0.5 + k) * DOS::step;
				this->fillHamiltonian(gamma);
				solver.compute(this->hamilton, false);

				energy += std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), double{},
					[this, k](double current, double toAdd) {
						return current + toAdd * BaseModel<DataType>::fermi_dirac(toAdd) * DOS::values[std::abs(k)];
					});
			}
			return 2 * energy / (2 * Constants::BASIS_SIZE - 1);
		};

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