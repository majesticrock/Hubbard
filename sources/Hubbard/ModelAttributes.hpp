#pragma once
#include "ModelParameters.hpp"
#include <complex>
#include <iostream>
#include <vector>
#include <algorithm>

namespace Hubbard {
	template <typename DataType>
	struct ModelAttributes {
	private:
		void initializeParamters(const ModelParameters& _params) {
			this->selfconsistency_values[0] = (std::abs(_params.U) + _params.V) * 0.5 + 0.1;
			this->selfconsistency_values[1] = std::abs(_params.U) * 0.5 + 0.1;
			this->selfconsistency_values[2] = std::abs(_params.U) * 0.5 + 0.1;
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
			this->selfconsistency_values[4] = std::abs(_params.V) * 0.2;

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
		ModelAttributes(ModelAttributes&& other) = default;
		ModelAttributes& operator=(const ModelAttributes& other) = default;
		ModelAttributes& operator=(ModelAttributes&& other) = default;

		explicit ModelAttributes(const ModelParameters& _params) {
			this->initializeParamters(_params);
		};

		ModelAttributes(const ModelAttributes<std::complex<double>>& other)
			: selfconsistency_values(other.selfconsistency_values.size()), converged{ other.converged }
		{
			if constexpr (std::is_floating_point_v<DataType>) {
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
				for (size_t i = 0U; i < selfconsistency_values.size(); ++i)
				{
					selfconsistency_values[i] = other.selfconsistency_values[i];
				}
			}
		};
		ModelAttributes(const ModelAttributes<double>& other)
			: selfconsistency_values(other.selfconsistency_values.size()), converged{ other.converged }
		{
			for (size_t i = 0U; i < selfconsistency_values.size(); ++i)
			{
				selfconsistency_values[i] = other.selfconsistency_values[i];
			}
		}

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
			for (const auto value : selfconsistency_values)
			{
				if (std::abs(value) > 1e-12) {
					return true;
				}
			}
			return false;
		};

		inline double renormalizedEnergy_up(const double GAMMA) const {
			if constexpr (std::is_same_v<DataType, std::complex<double>>) {
				return -(2. + this->selfconsistency_values[6].real()) * GAMMA;
			}
			else {
				return -(2. + this->selfconsistency_values[6]) * GAMMA;
			}
		};
		inline double renormalizedEnergy_down(const double GAMMA) const {
			if constexpr (std::is_same_v<DataType, std::complex<double>>) {
				return -(2. + this->selfconsistency_values[7].real()) * GAMMA;
			}
			else {
				return -(2. + this->selfconsistency_values[7]) * GAMMA;
			}
		};

		// Returns the total gap value sqrt(sc^2 + cdw^2)
		inline double getTotalGapValue() const {
			return sqrt(std::abs(selfconsistency_values[0]) * std::abs(selfconsistency_values[0])
				+ std::abs(selfconsistency_values[2]) * std::abs(selfconsistency_values[2]));
		};

		inline bool isFinite(size_t i, double epsilon = 1e-12) const {
			return (std::abs(selfconsistency_values[i]) > epsilon);
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
		inline ModelAttributes& operator*=(const double rhs) {
			for (size_t i = 0U; i < this->selfconsistency_values.size(); ++i)
			{
				this->selfconsistency_values[i] *= rhs;
			}
			return *this;
		};
		inline ModelAttributes& operator/=(const double rhs) {
			for (size_t i = 0U; i < this->selfconsistency_values.size(); ++i)
			{
				this->selfconsistency_values[i] /= rhs;
			}
			return *this;
		};
	};

	template <typename DataType>
	inline ModelAttributes<DataType> operator+(ModelAttributes<DataType> lhs, const ModelAttributes<DataType>& rhs) {
		return lhs += rhs;
	};

	template <typename DataType>
	inline ModelAttributes<DataType> operator*(ModelAttributes<DataType> lhs, double rhs) {
		return lhs *= rhs;
	};

	template <typename DataType>
	inline ModelAttributes<DataType> operator*(double lhs, ModelAttributes<DataType> rhs) {
		return rhs *= lhs;
	};

	template <typename DataType>
	inline ModelAttributes<DataType> operator/(ModelAttributes<DataType> lhs, double rhs) {
		return lhs /= rhs;
	};
}