#pragma once
#include "BaseModel.hpp"
#include "../../../FermionCommute/sources/Coefficient.hpp"
#include "NumericalMomentum.hpp"

namespace Hubbard {
	// maps an index; [0, N_K) -> [-pi, pi)
	template <typename T>
	inline double index_to_k_vector(const T index) {
		return (((index * L_PI) / Constants::K_DISCRETIZATION) - L_PI);
	};

	template <typename DataType, size_t Dimension>
	class MomentumBasedModel : public BaseModel<DataType>
	{
	protected:
		virtual void fillHamiltonian(const NumericalMomentum<Dimension>& k_values) = 0;
		virtual void addToParameterSet(const SpinorMatrix& rho, ComplexParameterVector& F, const NumericalMomentum<Dimension>& k_values) = 0;
	public:
		MomentumBasedModel(const ModelParameters& _params)
			: BaseModel<DataType>(_params) {};

		template<typename StartingValuesDataType>
		MomentumBasedModel(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: BaseModel<DataType>(_params, startingValues) {};

		inline double computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector2i& momentum) const {
			if (coeff.name == "\\epsilon_0") {
				NumericalMomentum<2> temp{index_to_k_vector(momentum(0)), index_to_k_vector(momentum(1))};
				return temp.unperturbed_energy() - this->chemical_potential;
			}
			if (coeff.name == "\\frac{U}{N}") {
				return this->U_OVER_N;
			}
			if (coeff.name == "\\tilde{V}") {
				NumericalMomentum<2> temp{index_to_k_vector(momentum(0)), index_to_k_vector(momentum(1))};
				return this->V_OVER_N * temp.gamma();
			}
			throw(std::invalid_argument("Could not find the coefficient: " + coeff.name));
		};
	};
}