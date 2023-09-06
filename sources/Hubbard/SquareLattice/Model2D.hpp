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

		virtual void computeExpectationValues(std::vector<MatrixCL>& expecs, std::vector<complex_prec>& sum_of_all) override {
			expecs = std::vector<MatrixCL>(8U, Matrix_L::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
			sum_of_all = std::vector<complex_prec>(8U, complex_prec{});

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					NumericalMomentum<2> ks{ k, l };
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