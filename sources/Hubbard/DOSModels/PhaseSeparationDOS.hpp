#include "DOSBasedModel.hpp"

namespace Hubbard{
	template <typename DataType, class DOS>
	class PhaseSeparationDOS : public DOSBasedModel<DataType, DOS>
	{
	private:
		int _extra_dimensions{};

		void init(){
			Constants::SPINOR_SIZE = 4 * _extra_dimensions;
			this->hamilton = SpinorMatrix::Zero(Constants::SPINOR_SIZE, Constants::SPINOR_SIZE);
			this->rho = SpinorMatrix::Zero(Constants::SPINOR_SIZE, Constants::SPINOR_SIZE);
		};

		inline void fillHamiltonian(const global_floating_type& gamma){
			DOSBasedModel<DataType, DOS>::fillHamiltonian(gamma);

			// For PS
			for(int extra = 0; extra < _extra_dimensions; ++extra){
				for (int i = 0; i < 4; ++i) {
					for (int j = 0; j < 4; ++j) {
						this->hamilton(4 * (1 + extra) + i, 4 * (1 + extra) + j) = this->hamilton(i, j);
					}
				}
				for (int i = 0; i < 4; ++i) {
					this->hamilton(4 * extra + i, 4 * (1 + extra) + i) = i < 2 ? DELTA_PS : -DELTA_PS;
					this->hamilton(4 * (1 + extra) + i, 4 * extra + i) = i < 2 ? DELTA_PS : -DELTA_PS;
				}
			}
		};

        inline void setParameterSet(ComplexParameterVector& F, const global_floating_type gamma) const {
			F.fill(complex_prec{});

            int[4] idx = {0, 1, 2, 3};
            for(int extra = 0; extra < _extra_dimensions; ++extra){
                for{int i = 0; i < 4; ++i}{
                    idx[i] = 4 * extra + i;
                }
                F(0) += -(this->rho(idx[0], idx[1]) + this->rho(idx[1], idx[0]) - this->rho(idx[2], idx[3]) - this->rho(idx[3], idx[2])).real(); // CDW
			    F(1) += -(this->rho(idx[0], idx[1]) + this->rho(idx[1], idx[0]) + this->rho(idx[2], idx[3]) + this->rho(idx[3], idx[2])).real(); // AFM
			    F(2) += -(this->rho(idx[0], idx[2]) + this->rho(idx[1], idx[3])); // SC
			    F(3) += -(2.0 / DOS::DIMENSION) * gamma * (this->rho(idx[0], idx[2]) - this->rho(idx[1], idx[3])); // Gamma SC
			    //F(4) += 0; // unused
			    F(5) += -(this->rho(idx[0], idx[3]) + this->rho(idx[1], idx[2])); // Eta
			    F(6) += -(2.0 / DOS::DIMENSION) * gamma * (this->rho(idx[0], idx[0]) - this->rho(idx[1], idx[1])).real(); // Gamma Occupation Up
			    F(7) += +(2.0 / DOS::DIMENSION) * gamma * (this->rho(idx[2], idx[2]) - this->rho(idx[3], idx[3])).real(); // Gamma Occupation Down

                if (extra + 1 < _extra_dimensions)
			        F(8) += -this->rho(idx[0], idx[0] + 4).real() - this->rho(idx[1], idx[1] + 4).real() + this->rho(idx[2], idx[2] + 4).real() + this->rho(idx[3], idx[3] + 4).real(); // PS
            }
            for(int i = 0; i < F.size() - 1; ++i){
                F(i) /= _extra_dimensions;
            }
            F(F.size() - 2) /= (_extra_dimensions - 1);
		};

        virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override {
			F.fill(global_floating_type{});
			std::conditional_t<Utility::is_complex<DataType>(), ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			std::copy(x.begin(), x.end(), this->model_attributes.begin());

			auto expectationValues = [this](global_floating_type gamma, ComplexParameterVector& result) {
				this->fillHamiltonian(gamma);
				this->fillRho();
				this->setParameterSet(result, gamma);
				};
	
			complex_F = this->_self_consistency_integrator.integrate_by_reference_symmetric(expectationValues);

			if constexpr (!std::is_same_v<DataType, complex_prec>) {
				complexParametersToReal(complex_F, F);
			}
			this->applyIteration(F);
			F -= x;
		};
	public:
		PhaseSeparationDOS(const ModelParameters& _params, int extra_dimensions) 
			: DOSBasedModel<DataType, DOS>(_params), _extra_dimensions(extra_dimensions)
		{
			init();
		};

		template<typename StartingValuesDataType>
		PhaseSeparationDOS(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues, int extra_dimensions)
			: DOSBasedModel<DataType, DOS>(_params, startingValues), _extra_dimensions(extra_dimensions)
		{
			init();
		};
	};
}