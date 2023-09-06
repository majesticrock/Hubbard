#pragma once
#include "BaseModel.hpp"
#include "../../../FermionCommute/sources/Coefficient.hpp"
#include "NumericalMomentum.hpp"

namespace Hubbard {
	template <typename DataType, size_t Dimension>
	class MomentumBasedModel : public BaseModel<DataType>
	{
	protected:
		using ParameterVector = typename BaseModel<DataType>::ParameterVector;
		virtual void fillHamiltonian(const NumericalMomentum<Dimension>& k_values) = 0;
		virtual void addToParameterSet(ComplexParameterVector& F, const NumericalMomentum<Dimension>& k_values) = 0;
	public:
		virtual void iterationStep(const ParameterVector& x, ParameterVector& F) override {
			F.fill(global_floating_type{});
			std::conditional_t<Utility::is_complex<DataType>(),
				ComplexParameterVector&, ComplexParameterVector> complex_F = F;

			std::copy(x.begin(), x.end(), this->model_attributes.begin());

			NumericalMomentum<Dimension> ks;
			while (ks.iterateHalfBZ()) {
				this->fillHamiltonian(ks);
				this->fillRho();
				this->addToParameterSet(complex_F, ks);
			}

			if constexpr (!Utility::is_complex<DataType>()) {
				complexParametersToReal(complex_F, F);
			}
			this->applyIteration(F);

			F -= x;
		};
		
		MomentumBasedModel(const ModelParameters& _params)
			: BaseModel<DataType>(_params) {};

		template<typename StartingValuesDataType>
		MomentumBasedModel(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: BaseModel<DataType>(_params, startingValues) {};

		inline global_floating_type computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector<int, Dimension>& momentum) const {
			if (coeff.name == "\\epsilon_0") {
				NumericalMomentum<2> temp{{momentum(0) - Constants::K_DISCRETIZATION, momentum(1) - Constants::K_DISCRETIZATION}};
				return temp.unperturbed_energy() - this->chemical_potential;
			}
			if (coeff.name == "\\frac{U}{N}") {
				return this->U_OVER_N;
			}
			if (coeff.name == "\\tilde{V}") {
				NumericalMomentum<2> temp{{momentum(0) - Constants::K_DISCRETIZATION, momentum(1) - Constants::K_DISCRETIZATION}};
				return this->V_OVER_N * temp.gamma();
			}
			throw(std::invalid_argument("Could not find the coefficient: " + coeff.name));
		};
	};
}