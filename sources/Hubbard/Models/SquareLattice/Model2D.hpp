#pragma once
#include "../MomentumBasedModel.hpp"

#define DELTA_CDW this->model_attributes[0]
#define DELTA_AFM this->model_attributes[1]
#define DELTA_SC this->model_attributes[2]
#define GAMMA_SC this->model_attributes[3]
#define XI_SC this->model_attributes[4]
#define DELTA_ETA this->model_attributes[5]
#define GAMMA_OCCSpinUpATION_SpinUp this->model_attributes[6]
#define GAMMA_OCCSpinUpATION_SpinDown this->model_attributes[7]

namespace Hubbard::Models::SquareLattice
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

		virtual void computeExpectationValues(std::vector<ValueArray>& expecs, ValueArray& sum_of_all) override {
			expecs = std::vector<ValueArray>(8U, ValueArray::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
			// The first column saves the plain sums sum_q <expec>. The second column contains 0.5 * sum_q gamma(q) <expec>
			sum_of_all = ValueArray::Zero(8, 2);

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					NumericalMomentum<2> ks{ k, l };
					this->fillHamiltonian(ks);
					this->fillRho();

					// n_up
					expecs[0](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = this->get_n_up();
					// g_up
					expecs[1](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = std::real(this->get_g_up());
					// f
					expecs[2](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = std::real(this->get_f());
					// eta
					expecs[3](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = std::real(this->get_eta());
					// n_down
					expecs[4](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = this->get_n_down();
					// g_down
					expecs[5](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = std::real(this->get_g_down());
					// n_up + n_down
					expecs[6](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = this->get_n_up_plus_down();
					// g_up + g_down
					expecs[7](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = std::real(this->get_g_up_plus_down());
					for (size_t idx = 0U; idx < 8U; ++idx)
					{
						sum_of_all(idx, 0) += expecs[idx](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION);
						sum_of_all(idx, 1) += 0.5 * ks.gamma() * expecs[idx](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION);
					}

					if (abs(this->rho(3, 0)) > 1e-10) {
						std::cerr << "Warning: <eta> does not vanish! " << this->rho(3, 0) << std::endl;
					}
				}
			}
		};
	};
}