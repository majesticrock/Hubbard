#pragma once
#include "ModelParameters.hpp"
#include "../Utility/IsComplex.hpp"
#include <complex>
#include <iostream>
#include <vector>
#include <algorithm>
#include <initializer_list>

namespace Hubbard {
	template <typename DataType>
	struct ModelAttributes {
	private:
		void initializeParamters_3d(const ModelParameters& _params) {
			auto guess = [&]() -> double {
				if (std::abs(_params.U) > 1e-12) {
					return std::abs(_params.U) * exp(0.5 * log(36.)) * exp(-2. / (0.288731210720569176L * std::abs(_params.U)));
				}
				return 0.0;
				};
			this->selfconsistency_values[0] = guess() + std::abs(_params.V);
			this->selfconsistency_values[1] = guess() + std::abs(_params.V);
			this->selfconsistency_values[2] = guess() + std::abs(_params.V);
			if (_params.U > 0) {
				this->selfconsistency_values[2] = 0.;
			}
			else {
				this->selfconsistency_values[1] = 0.;
			}
			if (_params.V > 0) {
				this->selfconsistency_values[2] = 0.;
			}
			else if (_params.V < 0) {
				this->selfconsistency_values[0] = 0.;
			}

			this->selfconsistency_values[3] = 0.;
			this->selfconsistency_values[4] = abs(_params.V) * 0.2;

			this->selfconsistency_values[5] = 0.;//_params.U * 0.1;
			this->selfconsistency_values[6] = _params.V * 0.2;
			this->selfconsistency_values[7] = _params.V * 0.2;
		};
		void initializeParamters_2d(const ModelParameters& _params) {
			auto guess = [&]() -> double {
				if (std::abs(_params.U) > 1e-12) {
					return std::abs(_params.U) * 4. * exp(-2 * 3.1415926L / sqrt(std::abs(_params.U)));
				}
				return 0.0;
				};
			this->selfconsistency_values[0] = guess() + std::abs(_params.V);
			this->selfconsistency_values[1] = guess() + std::abs(_params.V);
			this->selfconsistency_values[2] = guess() + std::abs(_params.V);
			if (_params.U >= 0) {
				this->selfconsistency_values[2] = 0.;
			}
			if (_params.U <= 0) {
				this->selfconsistency_values[1] = 0.;
			}
			if (_params.V > 0) {
				this->selfconsistency_values[2] = 0.;
			}
			else if (_params.V < 0) {
				this->selfconsistency_values[0] = 0.;
			}

			this->selfconsistency_values[3] = 0.;
			this->selfconsistency_values[4] = abs(_params.V) * 0.2;

			this->selfconsistency_values[5] = 0.;//_params.U * 0.1;
			this->selfconsistency_values[6] = _params.V * 0.2;
			this->selfconsistency_values[7] = _params.V * 0.2;
		};
		void initializeParamters(const ModelParameters& _params) {
			this->selfconsistency_values[0] = (abs(_params.U) + _params.V) * 0.5 + 0.1;
			this->selfconsistency_values[1] = abs(_params.U) * 0.5 + 0.1;
			this->selfconsistency_values[2] = abs(_params.U) * 0.5 + 0.1;
			if (_params.U >= 0) {
				this->selfconsistency_values[2] = 0.;
			}
			if (_params.U <= 0) {
				this->selfconsistency_values[1] = 0.;
			}
			if (_params.V > 0) {
				this->selfconsistency_values[2] = 0.;
			}
			else if (_params.V < 0) {
				this->selfconsistency_values[0] = 0.;
			}

			this->selfconsistency_values[3] = 0.;
			this->selfconsistency_values[4] = abs(_params.V) * 0.2;

			this->selfconsistency_values[5] = 0.;//_params.U * 0.1;
			this->selfconsistency_values[6] = _params.V * 0.2;
			this->selfconsistency_values[7] = _params.V * 0.2;
		}
	public:
		/* On a 2D System we choose the following:
		* 0 - Delta_CDW
		* 1 - Delta_AFM
		* 2 - Delta_SC
		* 3 - Gamma_SC
		* 4 - Xi_SC (1D: Tau_SC, 3D: None? P-wave?)
		* 5 - Delta_eta
		* 6 - Gamma_n_up
		* 7 - Gamma_n_down
		*/
		std::vector<DataType> selfconsistency_values = std::vector<DataType>(8);
		bool converged{};

		~ModelAttributes() = default;
		ModelAttributes() = default;
		ModelAttributes(std::initializer_list<DataType> i_list) : selfconsistency_values(i_list) {};
		ModelAttributes(ModelAttributes&& other) = default;
		ModelAttributes& operator=(const ModelAttributes& other) = default;
		ModelAttributes& operator=(ModelAttributes&& other) = default;

		explicit ModelAttributes(const ModelParameters& _params, int dimension = 0) {
			if (dimension == 2) {
				this->initializeParamters_2d(_params);
			}
			else if (dimension == 3) {
				this->initializeParamters_3d(_params);
			}
			else {
				this->initializeParamters(_params);
			}
		};

		ModelAttributes(const ModelAttributes<complex_prec>& other)
			: selfconsistency_values(other.selfconsistency_values.size()), converged{ other.converged }
		{
			if constexpr (!Utility::is_complex<DataType>()) {
				for (size_t i = 0U; i < selfconsistency_values.size(); ++i)
				{
					if (i == 4 || i == 7) {
						selfconsistency_values[i] = other.selfconsistency_values[i].imag();
					}
					else {
						selfconsistency_values[i] = other.selfconsistency_values[i].real();
					}
				}
			}
			else {
				std::copy(other.begin(), other.end(), this->begin());
			}
		};
		ModelAttributes(const ModelAttributes<global_floating_type>& other)
			: selfconsistency_values(other.selfconsistency_values.size()), converged{ other.converged }
		{
			std::copy(other.begin(), other.end(), this->begin());
		};

		/*
		* Utility functions
		*/

		inline DataType& operator[](size_t i) {
			assert(i < selfconsistency_values.size());
			return selfconsistency_values[i];
		};
		inline const DataType& operator[](size_t i) const {
			assert(i < selfconsistency_values.size());
			return selfconsistency_values[i];
		};
		inline size_t size() const noexcept {
			return selfconsistency_values.size();
		}
		inline void push_back(const DataType& value) {
			selfconsistency_values.push_back(value);
		};
		inline void push_back(DataType&& value) {
			selfconsistency_values.push_back(std::move(value));
		};
		inline auto begin() {
			return selfconsistency_values.begin();
		}
		inline auto begin() const {
			return selfconsistency_values.begin();
		}
		inline auto end() {
			return selfconsistency_values.end();
		}
		inline auto end() const {
			return selfconsistency_values.end();
		}
		inline void reset() {
			converged = false;
			std::fill(begin(), end(), DataType{});
		}

		inline bool isOrdered() const {
			for (const auto& value : selfconsistency_values)
			{
				if (abs(value) > 1e-12) {
					return true;
				}
			}
			return false;
		};

		inline global_floating_type renormalizedEnergy_up(const global_floating_type GAMMA) const {
			if constexpr (Utility::is_complex<DataType>()) {
				return -(2. + this->selfconsistency_values[6].real()) * GAMMA;
			}
			else {
				return -(2. + this->selfconsistency_values[6]) * GAMMA;
			}
		};
		inline global_floating_type renormalizedEnergy_down(const global_floating_type GAMMA) const {
			if constexpr (Utility::is_complex<DataType>()) {
				return -(2. + this->selfconsistency_values[7].real()) * GAMMA;
			}
			else {
				return -(2. + this->selfconsistency_values[7]) * GAMMA;
			}
		};

		// Returns the total gap value sqrt(sc^2 + cdw^2)
		inline global_floating_type getTotalGapValue() const {
			if constexpr (Utility::is_complex<decltype(selfconsistency_values[0])>()) {
				return sqrt(std::conj(selfconsistency_values[0]) * selfconsistency_values[0]
					+ std::conj(selfconsistency_values[1]) * selfconsistency_values[1]
					+ std::conj(selfconsistency_values[2]) * selfconsistency_values[2]);
			}
			else {
				return sqrt(selfconsistency_values[0] * selfconsistency_values[0]
					+ selfconsistency_values[1] * selfconsistency_values[1]
					+ selfconsistency_values[2] * selfconsistency_values[2]);
			}
		};

		inline bool isFinite(const size_t i, const double epsilon = 1e-12) const {
			return (abs(selfconsistency_values[i]) > epsilon);
		}
		inline void print() const {
			for (const auto& value : selfconsistency_values)
			{
				std::cout << value << "\t";
			}

			std::cout << "\n    Delta_tot = " << getTotalGapValue() << std::endl;
		};

		/*
		* Arithmetric operators
		*/

		inline ModelAttributes& operator+=(const ModelAttributes& rhs) {
			for (size_t i = 0U; i < this->selfconsistency_values.size(); ++i)
			{
				this->selfconsistency_values[i] += rhs.selfconsistency_values[i];
			}
			return *this;
		};
		template <class RealType>
		inline ModelAttributes& operator*=(const RealType rhs) {
			for (auto& value : this->selfconsistency_values)
			{
				value *= rhs;
			}
			return *this;
		};
		template <class RealType>
		inline ModelAttributes& operator/=(const RealType rhs) {
			for (auto& value : this->selfconsistency_values)
			{
				value /= rhs;
			}
			return *this;
		};
	};

	template <typename DataType>
	inline ModelAttributes<DataType> operator+(ModelAttributes<DataType> lhs, const ModelAttributes<DataType>& rhs) {
		return lhs += rhs;
	};

	template <typename DataType, class RealType>
	inline ModelAttributes<DataType> operator*(ModelAttributes<DataType> lhs, RealType rhs) {
		return lhs *= rhs;
	};

	template <typename DataType, class RealType>
	inline ModelAttributes<DataType> operator*(RealType lhs, ModelAttributes<DataType> rhs) {
		return rhs *= lhs;
	};

	template <typename DataType, class RealType>
	inline ModelAttributes<DataType> operator/(ModelAttributes<DataType> lhs, RealType rhs) {
		return lhs /= rhs;
	};

	template <typename DataType>
	inline std::ostream& operator<<(std::ostream& os, ModelAttributes<DataType> const& attributes) {
		for (const auto& value : attributes) {
			os << value << "  ";
		}
		return os;
	};
}