#pragma once
#include "BaseModel.hpp"
#include "../../../FermionCommute/sources/Coefficient.hpp"
#include <array>
#include <numeric>

namespace Hubbard {
	template<size_t Dimension>
	inline double tau(const std::array<double, Dimension>& ks){
		return std::accumulate(ks.begin(), ks.end(), double{}, [](double current, double toAdd){
			return current + sin(toAdd);
		});
	}
	template<size_t Dimension>
	inline double gamma(const std::array<double, Dimension>& ks){
		return std::accumulate(ks.begin(), ks.end(), double{}, [](double current, double toAdd){
			return current + cos(toAdd);
		});
	}
	inline double xi(double k_x, double k_y) {
		return cos(k_x) - cos(k_y);
	}
	inline double xi(const std::array<double, 2>& ks){
		return xi(ks[0], ks[1]);
	}
	template<size_t Dimension>
	inline double unperturbed_energy(const std::array<double, Dimension>& ks){
		return -2. * gamma(ks);
	}
	// maps an index; [0, N_K) -> [-pi, pi)
	template <typename T>
	inline double index_to_k_vector(const T index) {
		return (((index * L_PI) / Constants::K_DISCRETIZATION) - L_PI);
	};

	template <typename DataType, size_t Dimension>
	class MomentumBasedModel : public BaseModel<DataType>
	{
	protected:
		virtual void fillHamiltonian(const std::array<double, Dimension>& k_values) = 0;
		virtual void addToParameterSet(const SpinorMatrix& rho, ComplexParameterVector& F, const std::array<double, Dimension>& k_values) = 0;
	public:
		MomentumBasedModel(const ModelParameters& _params) 
			: BaseModel<DataType>(_params) {};

		template<typename StartingValuesDataType>
		MomentumBasedModel(const ModelParameters& _params, const ModelAttributes<StartingValuesDataType>& startingValues)
			: BaseModel<DataType>(_params, startingValues) {};

		inline double computeCoefficient(const SymbolicOperators::Coefficient& coeff, const std::array<double, Dimension>& momentum) const {
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