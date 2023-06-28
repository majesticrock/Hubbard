#pragma once
#include "BaseModel.hpp"
#include "../../../FermionCommute/sources/Coefficient.hpp"
#include "BaseModelAttributes.hpp"

#define UNPACK_MOMENTUM(name) double_prec name = va_arg(args, double_prec);
#define UNPACK_1D UNPACK_MOMENTUM(k_x)
#define UNPACK_2D UNPACK_1D UNPACK_MOMENTUM(k_y)
#define UNPACK_3D UNPACK_2D UNPACK_MOMENTUM(k_z)

namespace Hubbard {
	template <typename DataType, int Dimension>
	class MomentumBasedModel : public BaseModel<DataType>
	{
	protected:
		// maps an index; [0, N_K) -> [-pi, pi)
		template <typename T>
		inline double_prec index_to_k_vector(const T index) const {
			return (((index * L_PI) / Constants::K_DISCRETIZATION) - L_PI);
		};
	public:
		MomentumBasedModel(const ModelParameters& _params) : BaseModel<DataType>(_params) {};
		MomentumBasedModel(const ModelParameters& _params, const typename BaseModel<DataType>::BaseAttributes& startingValues)
			: BaseModel<DataType>(_params, startingValues) {};

		virtual inline double_prec computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector<int, Dimension>& momentum) const {
			if (coeff.name == "\\epsilon_0") {
				return (unperturbed_energy(index_to_k_vector(momentum(0)), index_to_k_vector(momentum(1))) - this->chemical_potential);
			}
			if (coeff.name == "\\frac{U}{N}") {
				return this->U_OVER_N;
			}
			if (coeff.name == "\\tilde{V}") {
				return this->V_OVER_N * (cos(index_to_k_vector(momentum(0))) + cos(index_to_k_vector(momentum(1))));
			}
			throw(std::invalid_argument("Could not find the coefficient: " + coeff.name));
		};
	};
}