#pragma once
#include "../MomentumBasedModel.hpp"

#define DELTA_CDW this->model_attributes[0]
#define DELTA_AFM this->model_attributes[1]
#define DELTA_SC this->model_attributes[2]
#define GAMMA_SC this->model_attributes[3]
#define XI_SC this->model_attributes[4]
#define DELTA_ETA this->model_attributes[5]
#define GAMMA_OCCUPATION_UP this->model_attributes[6]
#define GAMMA_OCCUPATION_DOWN this->model_attributes[7]

namespace Hubbard::SquareLattice
{
	inline global_floating_type xi(global_floating_type k_x, global_floating_type k_y) {
		return cos(k_x) - cos(k_y);
	}
	inline global_floating_type xi(const NumericalMomentum<2>& ks) {
		return xi(ks[0], ks[1]);
	}

	template <typename DataType>
	class Model2D : public MomentumBasedModel<DataType, 2>
	{
	protected:
		using ParameterVector = typename BaseModel<DataType>::ParameterVector;
	public:
		Model2D(const ModelParameters& _params) : MomentumBasedModel<DataType, 2>(_params) { };

		template<typename StartingValuesDataType>
		Model2D(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues) : MomentumBasedModel<DataType, 2>(_params, startingValues) { };

		// saves all one particle energies to reciever
		void getAllEnergies(std::vector<std::vector<global_floating_type>>& reciever)
		{
			reciever.resize(4 * Constants::K_DISCRETIZATION, std::vector<global_floating_type>(2 * Constants::K_DISCRETIZATION));
			Eigen::SelfAdjointEigenSolver<SpinorMatrix> solver;
			global_floating_type k_val{};
			global_floating_type l_val{};
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				k_val = k * BASE_PI / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					l_val = l * BASE_PI / Constants::K_DISCRETIZATION;
					NumericalMomentum<2> ks{k_val, l_val};
					this->fillHamiltonian(ks);
					solver.compute(this->hamilton, false);
					reciever[k + Constants::K_DISCRETIZATION][l + Constants::K_DISCRETIZATION] = solver.eigenvalues()(0);

					for (int i = 1; i < 4; i++)
					{
						if (abs(solver.eigenvalues()(0) - solver.eigenvalues()(i)) > 1e-8) {
							reciever[k + 3 * Constants::K_DISCRETIZATION][l + Constants::K_DISCRETIZATION] = solver.eigenvalues()(i);
							break;
						}
					}
				}
			}
		};

		inline virtual global_floating_type entropyPerSite() override {
			using std::log;
			global_floating_type entropy{};
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; ++k)
			{
				global_floating_type k_x = (k * BASE_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < 0; ++l)
				{
					global_floating_type k_y = (l * BASE_PI) / Constants::K_DISCRETIZATION;
					NumericalMomentum<2> ks{k_x, k_y};

					this->fillHamiltonian(ks);
					this->hamilton_solver.compute(this->hamilton, false);

					entropy += std::accumulate(this->hamilton_solver.eigenvalues().begin(), this->hamilton_solver.eigenvalues().end(), global_floating_type{},
						[this](global_floating_type current, global_floating_type toAdd) {
							auto occ = BaseModel<DataType>::fermi_dirac(toAdd);
							// Let's just not take the ln of 0. Negative numbers cannot be reached (because math...)
							return (occ > 1e-12 ? current - occ * log(occ) : current);
						});
				}
			}
			return entropy / Constants::BASIS_SIZE;
		};

		inline virtual global_floating_type internalEnergyPerSite() override {
			global_floating_type energy{};
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; ++k)
			{
				global_floating_type k_x = (k * BASE_PI) / Constants::K_DISCRETIZATION;
				for (int l = -Constants::K_DISCRETIZATION; l < 0; ++l)
				{
					global_floating_type k_y = (l * BASE_PI) / Constants::K_DISCRETIZATION;
					NumericalMomentum<2> ks{k_x, k_y};
					this->fillHamiltonian(ks);
					this->hamilton_solver.compute(this->hamilton, false);

					energy += std::accumulate(this->hamilton_solver.eigenvalues().begin(), this->hamilton_solver.eigenvalues().end(), global_floating_type{},
						[this](global_floating_type current, global_floating_type toAdd) {
							return current + toAdd * BaseModel<DataType>::fermi_dirac(toAdd);
						});
				}
			}
			return energy / Constants::BASIS_SIZE;
		};

		inline virtual void computeExpectationValues(std::vector<MatrixCL>& expecs, std::vector<complex_prec>& sum_of_all) override {
			expecs = std::vector<MatrixCL>(8U, Matrix_L::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
			sum_of_all = std::vector<complex_prec>(8U, complex_prec{});

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					NumericalMomentum<2> ks{ (k* BASE_PI) / Constants::K_DISCRETIZATION, (l* BASE_PI) / Constants::K_DISCRETIZATION };
					this->fillHamiltonian(ks);
					this->fillRho();

					// n_up
					expecs[0](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = 1 - this->rho(0, 0).real();
					// g_up
					expecs[1](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = -this->rho(1, 0);
					// f
					expecs[2](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = -this->rho(0, 2);
					// eta
					expecs[3](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = -this->rho(0, 3);
					// n_down
					expecs[4](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = this->rho(2, 2).real();
					// g_down
					expecs[5](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = this->rho(2, 3);
					// n_up + n_down
					expecs[6](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = 1 - (this->rho(0, 0) - this->rho(2, 2)).real();
					// g_up + g_down
					expecs[7](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = this->rho(2, 3) - this->rho(1, 0);
					for (size_t idx = 0U; idx < 8U; ++idx)
					{
						sum_of_all[idx] += expecs[idx](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION);
					}

					if (abs(this->rho(3, 0)) > 1e-10) {
						std::cerr << "Warning: <eta> does not vanish! " << this->rho(3, 0) << std::endl;
					}
				}
			}
		};
	};
}