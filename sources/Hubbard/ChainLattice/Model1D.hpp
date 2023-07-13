#pragma once
#include "../MomentumBasedModel.hpp"
#include <numeric>

#define DELTA_CDW this->model_attributes[0]
#define DELTA_AFM this->model_attributes[1]
#define DELTA_SC this->model_attributes[2]
#define GAMMA_SC this->model_attributes[3]
#define TAU_SC this->model_attributes[4]
#define DELTA_ETA this->model_attributes[5]
#define GAMMA_OCCUPATION_UP this->model_attributes[6]
#define GAMMA_OCCUPATION_DOWN this->model_attributes[7]

namespace Hubbard::ChainLattice
{
	inline void complexParametersToReal(const ComplexParameterVector& c, Eigen::VectorXd& r) {
		r(0) = c(0).real(); // CDW
		r(1) = c(1).real(); // AFM
		r(2) = c(2).real(); // SC
		r(3) = c(3).real(); // Gamma SC
		r(4) = c(4).imag(); // Xi SC
		r(5) = c(5).imag(); // Eta
		r(6) = c(6).real(); // Gamma Occupation Up
		r(7) = c(7).real(); // Gamma Occupation Down
	};

	template <typename DataType>
	class Model1D : public MomentumBasedModel<DataType, 1>
	{
	protected:
		using ParameterVector = typename BaseModel<DataType>::ParameterVector;

		inline void iterationStep(const ParameterVector& x, ParameterVector& F) {
			F.fill(0);
			std::conditional_t<std::is_same_v<DataType, complex_prec>,
				ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			SpinorMatrix rho{ SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE) };
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			std::copy(x.begin(), x.end(), this->model_attributes.selfconsistency_values.begin());

			for (int k = -Constants::K_DISCRETIZATION; k < 0; ++k)
			{
				const NumericalMomentum<1> k_x{ (k* L_PI) / Constants::K_DISCRETIZATION };
				this->fillHamiltonian(NumericalMomentum{ k_x });
				solver.compute(this->hamilton);

				this->fillRho(rho, solver);
				this->addToParameterSet(rho, complex_F, k_x);
			}

			if constexpr (!std::is_same_v<DataType, complex_prec>) {
				complexParametersToReal(complex_F, F);
			}
			this->setParameters(F);
			F -= x;
		};
	public:
		Model1D(const ModelParameters& _params) : MomentumBasedModel<DataType, 1>(_params) {};

		template<typename StartingValuesDataType>
		Model1D(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: MomentumBasedModel<DataType, 1>(_params, startingValues) {};

		inline virtual double entropyPerSite() override {
			double entropy = 0;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				const NumericalMomentum<1> k_x{ (k* L_PI) / Constants::K_DISCRETIZATION };

				this->fillHamiltonian(k_x);
				solver.compute(this->hamilton, false);

				entropy += std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), double{},
					[this](double current, double toAdd) {
						auto occ = BaseModel<DataType>::fermi_dirac(toAdd);
						// Let's just not take the ln of 0. Negative numbers cannot be reached (because math...)
						return (occ > 1e-12 ? current - occ * std::log(occ) : current);
					});
			}
			return entropy / Constants::BASIS_SIZE;
		};

		inline virtual double internalEnergyPerSite() override {
			double energy = 0;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				const NumericalMomentum<1> k_x{ (k* L_PI) / Constants::K_DISCRETIZATION };

				this->fillHamiltonian(k_x);
				solver.compute(this->hamilton, false);

				energy += std::accumulate(solver.eigenvalues().begin(), solver.eigenvalues().end(), double{},
					[this](double current, double toAdd) {
						return current + toAdd * BaseModel<DataType>::fermi_dirac(toAdd);
					});
			}
			return energy / Constants::BASIS_SIZE;
		};

		inline void computeExpectationValues(std::vector<VectorCL>& expecs, std::vector<complex_prec>& sum_of_all) {
			SpinorMatrix rho = SpinorMatrix::Zero(this->SPINOR_SIZE, this->SPINOR_SIZE);
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;

			expecs = std::vector<VectorCL>(8, VectorCL::Zero(2 * Constants::K_DISCRETIZATION));
			sum_of_all = std::vector<std::complex<double>>(8, 0.0);

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				const NumericalMomentum<1> k_x{ (k* L_PI) / Constants::K_DISCRETIZATION };
				this->fillHamiltonian(k_x);
				solver.compute(this->hamilton);
				this->fillRho(rho, solver);

				// n_up
				expecs[0](k + Constants::K_DISCRETIZATION) = 1 - rho(0, 0).real();
				// g_up
				expecs[1](k + Constants::K_DISCRETIZATION) = -rho(1, 0);
				// f
				expecs[2](k + Constants::K_DISCRETIZATION) = -rho(0, 2);
				// eta
				expecs[3](k + Constants::K_DISCRETIZATION) = -rho(0, 3);
				// n_down
				expecs[4](k + Constants::K_DISCRETIZATION) = rho(2, 2).real();
				// g_down
				expecs[5](k + Constants::K_DISCRETIZATION) = rho(2, 3);
				// n_up + n_down
				expecs[6](k + Constants::K_DISCRETIZATION) = 1 - (rho(0, 0) - rho(2, 2)).real();
				// g_up + g_down
				expecs[7](k + Constants::K_DISCRETIZATION) = rho(2, 3) - rho(1, 0);
				for (int idx = 0; idx < 8; idx++)
				{
					sum_of_all[idx] += expecs[idx](k + Constants::K_DISCRETIZATION);
				}

				if (std::abs(rho(3, 0)) > 1e-10) {
					std::cerr << "Warning: <eta> does not vanish! " << rho(3, 0) << std::endl;
				}
			}
		};
	};
}