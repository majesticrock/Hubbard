#pragma once
#include "../MomentumBasedModel.hpp"
#include <type_traits>

namespace Hubbard::ChainLattice {
	template <typename DataType>
	class Model1D : public MomentumBasedModel<DataType, 1>
	{
	protected:
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

			for (size_t i = 0; i < this->parameterMapper.size(); i++)
			{
				*(this->parameterMapper[i]) = x(i);
			}

			for (int k = -Constants::K_DISCRETIZATION; k < 0; k++)
			{
				double_prec k_x = (k * L_PI) / Constants::K_DISCRETIZATION;
				this->fillHamiltonian(1, k_x);
				solver.compute(this->hamilton);

				this->fillRho(rho, solver);
				this->addToParameterSet(rho, complex_F, 1, k_x);
			}

			if constexpr (!std::is_same<DataType, complex_prec>::value) {
				complexParametersToReal(complex_F, F);
			}
			this->setParameters(F);
			F -= x;
		};
	public:
		Model1D(const ModelParameters& _params) : MomentumBasedModel<DataType, 1>(_params) {};
		Model1D(const ModelParameters& _params, const typename BaseModel<DataType>::BaseAttributes& startingValues)
			: MomentumBasedModel<DataType, 1>(_params, startingValues) {};

		inline virtual double_prec entropyPerSite() override {
			double_prec entropy = 0;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * L_PI) / Constants::K_DISCRETIZATION;

				this->fillHamiltonian(1, k_x);
				solver.compute(this->hamilton, false);

				for (size_t i = 0; i < solver.eigenvalues().size(); i++)
				{
					auto occ = BaseModel<DataType>::fermi_dirac(solver.eigenvalues()(i));
					if (occ > 1e-12) { // Let's just not take the ln of 0. Negative numbers cannot be reached (by definition)
						entropy -= occ * std::log(occ);
					}
				}
			}
			return entropy / Constants::BASIS_SIZE;
		};

		inline virtual double_prec internalEnergyPerSite() override {
			double_prec energy = 0;
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				double_prec k_x = (k * L_PI) / Constants::K_DISCRETIZATION;

				this->fillHamiltonian(1, k_x);
				solver.compute(this->hamilton, false);

				for (size_t i = 0; i < solver.eigenvalues().size(); i++)
				{
					energy += BaseModel<DataType>::fermi_dirac(solver.eigenvalues()(i)) * solver.eigenvalues()(i);
				}
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
				this->fillHamiltonian(1, (k * L_PI) / Constants::K_DISCRETIZATION);
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