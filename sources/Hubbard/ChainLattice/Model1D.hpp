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
	inline void complexParametersToReal(const ComplexParameterVector& c, Vector_L& r) {
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

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override {
			F.fill(global_floating_type{});
			std::conditional_t<Utility::is_complex<DataType>(),
				ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			std::copy(x.begin(), x.end(), this->model_attributes.begin());

			for (int k = -Constants::K_DISCRETIZATION; k < 0; ++k)
			{
				const NumericalMomentum<1> k_x{ k };

				this->fillHamiltonian(NumericalMomentum{ k_x });
				this->fillRho();
				this->addToParameterSet(complex_F, k_x);
			}

			if constexpr (!Utility::is_complex<DataType>()) {
				complexParametersToReal(complex_F, F);
			}
			this->applyIteration(F);
			F -= x;
		};
	public:
		Model1D(const ModelParameters& _params) : MomentumBasedModel<DataType, 1>(_params) {};

		template<typename StartingValuesDataType>
		Model1D(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: MomentumBasedModel<DataType, 1>(_params, startingValues) {};

		virtual void computeExpectationValues(std::vector<MatrixCL>& expecs, std::vector<complex_prec>& sum_of_all) override {
			expecs = std::vector<MatrixCL>(8, Matrix_L::Zero(2 * Constants::K_DISCRETIZATION, 1));
			sum_of_all = std::vector<complex_prec>(8, complex_prec{});

			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				const NumericalMomentum<1> k_x{ k };
				this->fillHamiltonian(k_x);
				this->fillRho();

				// n_up
				expecs[0](k + Constants::K_DISCRETIZATION, 0) = 1 - this->rho(0, 0).real();
				// g_up
				expecs[1](k + Constants::K_DISCRETIZATION, 0) = -this->rho(1, 0);
				// f
				expecs[2](k + Constants::K_DISCRETIZATION, 0) = -this->rho(0, 2);
				// eta
				expecs[3](k + Constants::K_DISCRETIZATION, 0) = -this->rho(0, 3);
				// n_down
				expecs[4](k + Constants::K_DISCRETIZATION, 0) = this->rho(2, 2).real();
				// g_down
				expecs[5](k + Constants::K_DISCRETIZATION, 0) = this->rho(2, 3);
				// n_up + n_down
				expecs[6](k + Constants::K_DISCRETIZATION, 0) = 1 - (this->rho(0, 0) - this->rho(2, 2)).real();
				// g_up + g_down
				expecs[7](k + Constants::K_DISCRETIZATION, 0) = this->rho(2, 3) - this->rho(1, 0);
				for (int idx = 0; idx < 8; idx++)
				{
					sum_of_all[idx] += expecs[idx](k + Constants::K_DISCRETIZATION, 0);
				}

				if (abs(this->rho(3, 0)) > 1e-10) {
					std::cerr << "Warning: <eta> does not vanish! " << this->rho(3, 0) << std::endl;
				}
			}
		};
	};
}