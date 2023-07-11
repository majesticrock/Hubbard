#pragma once
#include "BaseModel.hpp"
#include "../../../FermionCommute/sources/Coefficient.hpp"

#define UNPACK_MOMENTUM(name) double name = va_arg(args, double);
#define UNPACK_1D UNPACK_MOMENTUM(k_x)
#define UNPACK_2D UNPACK_1D UNPACK_MOMENTUM(k_y)
#define UNPACK_3D UNPACK_2D UNPACK_MOMENTUM(k_z)

namespace Hubbard {
	template<typename... Args>
	inline double gamma(Args... ks) {
		return (cos(ks) + ...);
	}
	inline double xi(double k_x, double k_y) {
		return cos(k_x) - cos(k_y);
	}
	template<typename... Args>
	inline double unperturbed_energy(Args... ks) {
		return -2. * (cos(ks) + ...);
	};
	// maps an index; [0, N_K) -> [-pi, pi)
	template <typename T>
	inline double index_to_k_vector(const T index) {
		return (((index * L_PI) / Constants::K_DISCRETIZATION) - L_PI);
	};

	template <typename DataType, int Dimension>
	class MomentumBasedModel : public BaseModel<DataType>
	{
	public:
		MomentumBasedModel(const ModelParameters& _params) 
			: BaseModel<DataType>(_params) {};

		template<typename StartingValuesDataType>
		MomentumBasedModel(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: BaseModel<DataType>(_params, startingValues) {};

		inline double computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector<int, Dimension>& momentum) const {
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