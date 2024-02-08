#include "DOSBasedModel.hpp"
#include "../Selfconsistency/BroydenSolver.hpp"

namespace Hubbard::DOSModels {
	template <class DOS>
	class PhaseSeparationDOS : public DOSBasedModel<global_floating_type, DOS>
	{
	private:
		using ParameterVector = typename BaseModel<global_floating_type>::ParameterVector;
		const size_t _MaxPreBroydenIterations;
		int _extra_dimensions{};

		void init() {
			Constants::SPINOR_SIZE = 4 + 4 * _extra_dimensions;
			this->hamilton = SpinorMatrix::Zero(Constants::SPINOR_SIZE, Constants::SPINOR_SIZE);
			this->rho = SpinorMatrix::Zero(Constants::SPINOR_SIZE, Constants::SPINOR_SIZE);
		};

		inline void fillHamiltonian(const global_floating_type& gamma) {
			DOSBasedModel<global_floating_type, DOS>::fillHamiltonian(gamma);

			// For PS
			for (int extra = 0; extra < _extra_dimensions; ++extra) {
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
			int idx[4] = { 0, 1, 2, 3 };

			F(0) += -(this->rho(idx[0], idx[1]) + this->rho(idx[1], idx[0]) - this->rho(idx[2], idx[3]) - this->rho(idx[3], idx[2])).real(); // CDW
			F(1) += -(this->rho(idx[0], idx[1]) + this->rho(idx[1], idx[0]) + this->rho(idx[2], idx[3]) + this->rho(idx[3], idx[2])).real(); // AFM
			F(2) += -(this->rho(idx[0], idx[2]) + this->rho(idx[1], idx[3])); // SC
			F(3) += -(2.0 / DOS::DIMENSION) * gamma * (this->rho(idx[0], idx[2]) - this->rho(idx[1], idx[3])); // Gamma SC
			//F(4) += 0; // unused
			F(5) += -(this->rho(idx[0], idx[3]) + this->rho(idx[1], idx[2])); // Eta
			F(6) += -(2.0 / DOS::DIMENSION) * gamma * (this->rho(idx[0], idx[0]) - this->rho(idx[1], idx[1])).real(); // Gamma Occupation Up
			F(7) += +(2.0 / DOS::DIMENSION) * gamma * (this->rho(idx[2], idx[2]) - this->rho(idx[3], idx[3])).real(); // Gamma Occupation Down

			F(8) += -this->rho(idx[0], idx[0] + 4).real() - this->rho(idx[1], idx[1] + 4).real() + this->rho(idx[2], idx[2] + 4).real() + this->rho(idx[3], idx[3] + 4).real(); // PS
		};

		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override {
			F.fill(global_floating_type{});
			ComplexParameterVector complex_F = F;

			std::copy(x.begin(), x.end(), this->model_attributes.begin());

			auto expectationValues = [this](global_floating_type gamma, ComplexParameterVector& result) {
				this->fillHamiltonian(gamma);
				this->fillRho();
				this->setParameterSet(result, gamma);
				};

			complex_F = this->_self_consistency_integrator.integrate_by_reference_symmetric(expectationValues);

			complexParametersToReal(complex_F, F);
			this->applyIteration(F);
			F -= x;
		};
	public:
		PhaseSeparationDOS(const ModelParameters& _params, int extra_dimensions, size_t MaxPreBroydenIterations = 300U)
			: DOSBasedModel<global_floating_type, DOS>(_params), _MaxPreBroydenIterations(MaxPreBroydenIterations), _extra_dimensions(extra_dimensions)
		{
			init();
		};

		template<typename StartingValuesglobal_floating_type>
		PhaseSeparationDOS(const ModelParameters& _params, const ModelAttributes<StartingValuesglobal_floating_type>& startingValues, int extra_dimensions, size_t MaxPreBroydenIterations = 300U)
			: DOSBasedModel<global_floating_type, DOS>(_params, startingValues), _MaxPreBroydenIterations(MaxPreBroydenIterations), _extra_dimensions(extra_dimensions)
		{
			init();
		};

		virtual ModelAttributes<global_floating_type> computePhases(const PhaseDebuggingPolicy debugPolicy = NoWarning) override
		{
			Selfconsistency::BroydenSolver solver(this, &this->model_attributes, _MaxPreBroydenIterations);
			return solver.computePhases(debugPolicy);
		};
	};
}