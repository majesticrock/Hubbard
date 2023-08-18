#pragma once
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <algorithm>

using boost::math::constants::pi;
using boost::math::constants::half_pi;
using boost::multiprecision::exp;
using boost::multiprecision::sinh;
using boost::multiprecision::cosh;
using boost::multiprecision::tanh;
using boost::multiprecision::pow;

namespace Hubbard::DensityOfStates {
	template <class T>
	void expand_vector(std::vector<T>& vec) {
		if (vec.size() == 0) {
			vec.push_back(T{});
			return;
		}
		int n = vec.size();
		vec.reserve(2 * n - 1);
		for (int i = n - 1; i > 0; --i)
		{
			vec.insert(vec.begin() + i, T{});
		}
	};

	template <class RealType, class FinalType = double>
	struct tanh_sinh_helper {
	private:
		using WeightType = FinalType;
		RealType _sinh_x{};

		// Standard from -1 to 1
		const FinalType _lower_border{ -1 };
		const FinalType _upper_border{ 1 };
		const FinalType _center{ (_upper_border + _lower_border) / 2 }; // (a+b)/2
		const FinalType _half_distance{ (_upper_border - _lower_border) / 2 }; // (b-a)/2

		FinalType _step{ 1 };
		size_t _level{};
		int _min_k{};

	public:
		struct SaveTo {
			std::vector<RealType>* abscissa;
			std::vector<RealType>* border_to_abscissa;
			std::vector<WeightType>* weights;
			std::vector<FinalType>* values;
		};

		inline void set_sinh_x(RealType x) {
			_sinh_x = sinh(x);
		};
		inline RealType compute_abscissa() const {
			return _upper_border - compute_upper_border_to_abscissa();// _center + _half_distance * std::tanh(LONG_PI_2 * _sinh_x);
		};
		// Computes b - x
		inline RealType compute_upper_border_to_abscissa() const {
			return _half_distance * 2. /
				(1 + exp(pi<RealType>() * _sinh_x));
		};
		// Computes a - x
		inline RealType compute_lower_border_to_abscissa() const {
			return -2. * _half_distance * exp(pi<RealType>() * _sinh_x) / (1 + exp(pi<RealType>() * _sinh_x));
		};
		inline WeightType compute_weight(int k) const {
			return (_step * half_pi<WeightType>() * cosh(k * _step))
				/ static_cast<WeightType>(pow(cosh(half_pi<WeightType>() * _sinh_x), 2));
		};

		inline FinalType half_distance() const {
			return _half_distance;
		};
		inline void increase_level(SaveTo& save_to) {
			++_level;
			_step /= 2;
			_min_k *= 2;

			expand_vector(*(save_to.abscissa));
			expand_vector(*(save_to.border_to_abscissa));
			expand_vector(*(save_to.weights));
			expand_vector(*(save_to.values));
		};
		inline FinalType step() const {
			return _step;
		};
		inline size_t level() const {
			return _level;
		};

		template< int cut_off, class UnaryFunction>
		FinalType initial_filling(const UnaryFunction& dos, SaveTo& save_to) {
			FinalType integral{};
			auto fill_vectors = [&](int k, bool sign) {
				do {
					set_sinh_x((sign ? k : -k) * _step);
					save_to.border_to_abscissa->push_back(compute_upper_border_to_abscissa());
					save_to.abscissa->push_back(compute_abscissa());
					save_to.weights->push_back(compute_weight((sign ? k : -k)));
					save_to.values->push_back(dos(save_to.abscissa->back(), save_to.border_to_abscissa->back()));
					++k;

					integral += save_to.values->back() * static_cast<FinalType>(save_to.weights->back());
				} while (abs(save_to.values->back() * static_cast<FinalType>(save_to.weights->back())) > boost::math::pow<cut_off>(10)
					|| abs(save_to.weights->back()) > 1e-8);
				return k - 1;
			};

			_min_k = -fill_vectors(0, false);

			std::reverse(save_to.values->begin(), save_to.values->end());
			std::reverse(save_to.weights->begin(), save_to.weights->end());
			std::reverse(save_to.abscissa->begin(), save_to.abscissa->end());
			std::reverse(save_to.border_to_abscissa->begin(), save_to.border_to_abscissa->end());

			fill_vectors(1, true);
			return integral * _half_distance;
		};

		template<class UnaryFunction>
		inline void compute_step(const UnaryFunction& dos, int k, SaveTo& save_to) {
			set_sinh_x((k + _min_k) * _step);
			save_to.abscissa->at(k) = compute_abscissa();
			save_to.border_to_abscissa->at(k) = compute_upper_border_to_abscissa();
			save_to.weights->at(k) = compute_weight(k + _min_k);
			save_to.values->at(k) = dos(save_to.abscissa->at(k), save_to.border_to_abscissa->at(k));
		};

		tanh_sinh_helper(RealType a = -1, RealType b = 1)
			: _lower_border(a), _upper_border(b), _center((b + a) / 2), _half_distance((b - a) / 2) {};
	};
}