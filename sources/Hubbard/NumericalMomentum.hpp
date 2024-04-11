#pragma once
#include <iterator>
#include <numeric>
#include <algorithm>
#include <array>
#include "GlobalDefinitions.hpp"
#include "Constants.hpp"

namespace Hubbard {
	inline int mod_two_pi(int index) {
		if (index >= -Constants::K_DISCRETIZATION && index < Constants::K_DISCRETIZATION)
			return index;
		if (index < 0)
			return ((1 + abs(index / TWO_K_DISC)) * TWO_K_DISC + index) % TWO_K_DISC - Constants::K_DISCRETIZATION;
		else
			return index % TWO_K_DISC - Constants::K_DISCRETIZATION;
	};

	template<unsigned int Dimension>
	struct NumericalMomentum {
		global_floating_type momenta[Dimension];
		int k[Dimension];

		NumericalMomentum()
		{
			this->reset();
		};
		explicit NumericalMomentum(int default_k) {
			std::fill(std::begin(momenta), std::end(momenta), Constants::PI_DIV_DISCRETIZATION * default_k);
			std::fill(std::begin(k), std::end(k), default_k);
		}
		template <typename... Args>
		explicit NumericalMomentum(Args... args) : k{ static_cast<int>(args)... } {
			static_assert(sizeof...(Args) == Dimension, "Incorrect number of arguments");
			static_assert(std::conjunction_v<std::is_integral<Args>...>, "All arguments must be integers");
			//((std::cout << ',' << std::forward<Args>(args)), ...);
			for (size_t d = 0U; d < Dimension; ++d)
			{
				momenta[d] = k[d] * Constants::PI_DIV_DISCRETIZATION;
			}
		}
		explicit NumericalMomentum(const Eigen::Array<int, Dimension, 1>& point_in_bz) {
			for (size_t d = 0U; d < Dimension; ++d)
			{
				k[d] = point_in_bz(d);
				momenta[d] = k[d] * Constants::PI_DIV_DISCRETIZATION;
			}
		};
		explicit NumericalMomentum(const Eigen::Vector<global_floating_type, Dimension>& momentum) {
			for (size_t d = 0U; d < Dimension; ++d)
			{
				momenta[d] = momentum(d);
			}
		};

		inline static NumericalMomentum<Dimension> GammaPoint() {
			return NumericalMomentum<Dimension>(0);
		};
		inline static NumericalMomentum<Dimension> Q() {
			return NumericalMomentum<Dimension>(-Constants::K_DISCRETIZATION);
		}

		inline global_floating_type tau() const {
			return std::accumulate(std::begin(momenta), std::end(momenta), global_floating_type{}, [](global_floating_type current, global_floating_type toAdd) {
				return current + sin(toAdd);
				});
		};
		inline global_floating_type gamma() const {
			return std::accumulate(std::begin(momenta), std::end(momenta), global_floating_type{}, [](global_floating_type current, global_floating_type toAdd) {
				return current + cos(toAdd);
				});
		};
		inline global_floating_type unperturbed_energy() const {
			return -2. * gamma();
		};

		inline global_floating_type& operator[](unsigned int index) {
			assert(index < Dimension);
			return momenta[index];
		};
		inline const global_floating_type& operator[](unsigned int index) const {
			assert(index < Dimension);
			return momenta[index];
		};

		inline NumericalMomentum& operator++() {
			_increment<0>();
			return *this;
		};
		// To be used with a do-while(.iterateHalfBZ()) loop
		// Otherwise the first point is omitted!
		inline bool iterateHalfBZ() {
			_increment<0>();
			return (k[Dimension - 1] < 0);
		};
		// To be used with a do-while(.iterateFullBZ()) loop
		// Otherwise the first point is omitted!
		inline bool iterateFullBZ() {
			_increment<0>();
			return (k[Dimension - 1] < Constants::K_DISCRETIZATION);
		};

		template <unsigned int d>
		inline bool increase_d() {
			static_assert(d < Dimension, "NumericalMomentum out of bounds!");
			momenta[d] = (++k[d]) * Constants::PI_DIV_DISCRETIZATION;
			return (k[d] < Constants::K_DISCRETIZATION);
		}

		inline size_t getIndex() const {
			assert(k[0] + Constants::K_DISCRETIZATION >= 0);
			size_t index{ static_cast<size_t>(k[0] + Constants::K_DISCRETIZATION) };
			for (unsigned int i = 1U; i < Dimension; ++i)
			{
				index += 2 * i * Constants::K_DISCRETIZATION * (k[i] + Constants::K_DISCRETIZATION);
			}
			return index;
		};
		inline void reset() {
			std::fill(std::begin(momenta), std::end(momenta), -BASE_PI);
			std::fill(std::begin(k), std::end(k), -Constants::K_DISCRETIZATION);
		};

		bool isZero() const {
			return std::all_of(std::begin(k), std::end(k), [](int value) { return value == 0; });
		};

		bool isQ() const {
			return std::all_of(std::begin(k), std::end(k), [](int value) { return value == -Constants::K_DISCRETIZATION; });
		};

		global_floating_type squared_norm() const {
			return std::accumulate(std::begin(momenta), std::end(momenta), global_floating_type{}, [](global_floating_type current, global_floating_type toAdd) {
				return current + toAdd * toAdd;
				});
		};

		global_floating_type norm() const {
			return sqrt(squared_norm());
		};

		// This function assumes that no anomalous momenta (i.e. outside of [-pi, pi)) are being used
		inline NumericalMomentum& operator+=(const NumericalMomentum& rhs) {
			for (size_t i = 0U; i < Dimension; ++i)
			{
				this->k[i] += rhs.k[i];
				if (this->k[i] < -Constants::K_DISCRETIZATION) this->k[i] += TWO_K_DISC;
				else if (this->k[i] >= Constants::K_DISCRETIZATION) this->k[i] -= TWO_K_DISC;

				this->momenta[i] = this->k[i] * Constants::PI_DIV_DISCRETIZATION;
			}
			return *this;
		};
		// This function assumes that no anomalous momenta (i.e. outside of [-pi, pi)) are being used
		inline NumericalMomentum& operator-=(const NumericalMomentum& rhs) {
			for (size_t i = 0U; i < Dimension; ++i)
			{
				this->k[i] -= rhs.k[i];
				if (this->k[i] < -Constants::K_DISCRETIZATION) this->k[i] += TWO_K_DISC;
				else if (this->k[i] >= Constants::K_DISCRETIZATION) this->k[i] -= TWO_K_DISC;

				this->momenta[i] = this->k[i] * Constants::PI_DIV_DISCRETIZATION;
			}
			return *this;
		};

		inline NumericalMomentum& operator*=(const int rhs) {
			for (size_t i = 0U; i < Dimension; ++i)
			{
				this->k[i] = mod_two_pi(rhs * this->k[i]);
				this->momenta[i] = this->k[i] * Constants::PI_DIV_DISCRETIZATION;
			}
			return *this;
		};

	private:
		template <int _d>
		inline void _increment() {
			static_assert(_d < Dimension, "Call to increment in NumericalMomentum provides a too high dimension.");
			if (++k[_d] >= Constants::K_DISCRETIZATION) {
				if constexpr (_d + 1 < Dimension) k[_d] = -Constants::K_DISCRETIZATION;
				if constexpr (_d + 1 < Dimension) {
					_increment<_d + 1>();
				}
			}
			momenta[_d] = k[_d] * Constants::PI_DIV_DISCRETIZATION;
		};
	};

	template <unsigned int Dimension>
	inline NumericalMomentum<Dimension> operator+(NumericalMomentum<Dimension> lhs, const NumericalMomentum<Dimension>& rhs) {
		return lhs += rhs;
	};
	template <unsigned int Dimension>
	inline NumericalMomentum<Dimension> operator-(NumericalMomentum<Dimension> lhs, const NumericalMomentum<Dimension>& rhs) {
		return lhs -= rhs;
	};
	template <unsigned int Dimension>
	inline NumericalMomentum<Dimension> operator*(NumericalMomentum<Dimension> lhs, const int rhs) {
		return lhs *= rhs;
	};
	template <unsigned int Dimension>
	inline NumericalMomentum<Dimension> operator*(const int lhs, NumericalMomentum<Dimension> rhs) {
		return rhs *= lhs;
	};
	template <unsigned int Dimension>
	inline NumericalMomentum<Dimension> operator-(NumericalMomentum<Dimension> rhs) {
		return rhs *= -1;
	};

	template <unsigned int Dimension>
	inline std::ostream& operator<<(std::ostream& os, const NumericalMomentum<Dimension>& momentum) {
		os << momentum.momenta[0];
		for (size_t d = 1U; d < Dimension; ++d)
		{
			os << " " << momentum.momenta[d];
		}
		return os;
	};
}