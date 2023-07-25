#pragma once
#include "../BaseModel.hpp"
#include "../../../../FermionCommute/sources/Coefficient.hpp"
#include "../DensityOfStates/BaseDOS.hpp"
#include "../DensityOfStates/Square.hpp"
#include <algorithm>
#include <mutex>

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
		static constexpr bool DOS_IS_SQUARE = std::is_same_v<DOS, DensityOfStates::Square>;
		// The square lattice DOS cannot include gamma = 0, due to the singularity
		static int gammaLoopUpperBoundary;
		static std::mutex dos_mutex;

		void init() {
			this->model_attributes[4] = 0.;

			if (!DOS::computed) {
				std::lock_guard<std::mutex> guard(dos_mutex);
				// Might have been changed by another thread
				if (!DOS::computed) {
					gammaLoopUpperBoundary = DOSBasedModel<DataType, DOS>::DOS_IS_SQUARE ? Constants::BASIS_SIZE : Constants::BASIS_SIZE + 1;
					DOS dos;
					dos.computeValues();
					std::cout << "1 - DOS-Norm = " << std::scientific << 1. - dos.getNorm() << std::endl;
				}
			}
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

		void addToParameterSet(ComplexParameterVector& F, const double gamma, const double dos_value) {
			F(0) -= (this->rho(0, 1) + this->rho(1, 0) - this->rho(2, 3) - this->rho(3, 2)).real() * dos_value; // CDW
			F(1) -= (this->rho(0, 1) + this->rho(1, 0) + this->rho(2, 3) + this->rho(3, 2)).real() * dos_value; // AFM
			F(2) -= (this->rho(0, 2) + this->rho(1, 3)) * dos_value; // SC
			F(3) -= gamma * (this->rho(0, 2) - this->rho(1, 3)) * dos_value; // Gamma SC

			F(5) -= (this->rho(0, 3) + this->rho(1, 2)) * dos_value; // Eta
			F(6) -= gamma * (this->rho(0, 0) - this->rho(1, 1)).real() * dos_value; // Gamma Occupation Up
			F(7) += gamma * (this->rho(2, 2) - this->rho(3, 3)).real() * dos_value; // Gamma Occupation Down
		};

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override {
			F.fill(0);
			std::conditional_t<std::is_same_v<DataType, complex_prec>,
				ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			std::copy(x.begin(), x.end(), this->model_attributes.begin());

			for (int k = -Constants::BASIS_SIZE; k < gammaLoopUpperBoundary; ++k)
			{
				double gamma = k * DOS::step;
				if constexpr (this->DOS_IS_SQUARE) {
					gamma += 0.5 * DOS::step;
				}
				this->fillHamiltonian(gamma);
				this->fillRho();
				this->addToParameterSet(complex_F, gamma, DOS::values[k + Constants::BASIS_SIZE]);
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
			for (int k = -Constants::BASIS_SIZE; k < gammaLoopUpperBoundary; ++k)
			{
				double gamma = k * DOS::step;
				if constexpr (this->DOS_IS_SQUARE) {
					gamma += 0.5 * DOS::step;
				}
				this->fillHamiltonian(gamma);
				solver.compute(this->hamilton, false);
				entropy += std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), double{},
					[this, k](double current, double toAdd) {
						auto occ = BaseModel<DataType>::fermi_dirac(toAdd);
						// Let's just not take the ln of 0. Negative numbers cannot be reached (because math...)
						return (occ > 1e-12 ? current - occ * std::log(occ) : current) * DOS::values[k + Constants::BASIS_SIZE];
					});
			}
			return entropy / Constants::BASIS_SIZE;
		};

		inline virtual double internalEnergyPerSite() override {
			double energy = 0;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			for (int k = -Constants::BASIS_SIZE; k < gammaLoopUpperBoundary; ++k)
			{
				double gamma = k * DOS::step;
				if constexpr (this->DOS_IS_SQUARE) {
					gamma += 0.5 * DOS::step;
				}
				this->fillHamiltonian(gamma);
				solver.compute(this->hamilton, false);
				energy += std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), double{},
					[this, k](double current, double toAdd) {
						return current + toAdd * BaseModel<DataType>::fermi_dirac(toAdd) * DOS::values[k + Constants::BASIS_SIZE];
					});
			}
			return energy / Constants::BASIS_SIZE;
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

	template <typename DataType, class DOS>
	std::mutex DOSBasedModel<DataType, DOS>::dos_mutex;

	template <typename DataType, class DOS>
	int DOSBasedModel<DataType, DOS>::gammaLoopUpperBoundary = 0;
}