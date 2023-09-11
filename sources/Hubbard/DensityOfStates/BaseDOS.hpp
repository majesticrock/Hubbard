#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "../GlobalDefinitions.hpp"
#include "DOSIntegrator.hpp"

namespace Hubbard::DensityOfStates {
	using dos_precision = long_double_t;
	typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<100>> abscissa_t;

	_CONST_LONG_FLOATING R_AT_2 = (LONG_PI - 4) / 8;
	_CONST_LONG_FLOATING R_AT_2_1 = (5 * LONG_PI / 64 - 0.25L);
#ifdef _BOOST_PRECISION
	_CONST_LONG_FLOATING CUT_OFF = 1e-20;
#else
	_CONST_LONG_FLOATING CUT_OFF = 1e-12;
#endif

	// Computes sqrt(1 - x^2)
	template <class RealType>
	inline RealType sqrt_1_minus_x_squared(RealType x) {
		return sqrt((1 - x) * (1 + x));
	}
	// Use this overload for better results, if |x| is close to 1
	// However, you need to provide an accurate reprensantion of 1+/-x yourself.
	template <class RealType>
	inline RealType sqrt_1_minus_x_squared(RealType x, RealType one_minus_x) {
		return (x < 0.25 ? sqrt(1 - x * x) : sqrt(one_minus_x * (1 + x)));
	}

	template <bool includesZero, class DataType>
	inline void symmetrizeVector(std::vector<DataType>& v) {
		std::reverse(v.begin(), v.end());
		if constexpr (includesZero) {
			v.insert(v.end(), v.rbegin() + 1, v.rend());
		}
		else {
			v.insert(v.end(), v.rbegin(), v.rend());
		}
	}

	template <class DOS>
	inline dos_precision computeNorm() {
		
		return typename DOS::Integrator<dos_precision>().integrate_by_value([](dos_precision) -> dos_precision { return 1.L; });
	}

	struct BaseDOS {
		// Contains the values for gamma in (-2, 0) as the the positive part is symmetric.
		// Open boundaries, as the dos at 0 contains a singularity for d=2
		static std::vector<dos_precision> values;
		static std::vector<abscissa_t> abscissa;
		static dos_precision step;
		static bool computed;
		static dos_precision integrateValues();
		static void printValues();

		virtual void computeValues() = 0;
		virtual ~BaseDOS() = default;

		inline static global_floating_type abscissa_v(size_t index) {
			return static_cast<global_floating_type>(index < abscissa.size() ? abscissa[index] : -abscissa[index - abscissa.size()]);
		};
		inline static size_t size() noexcept {
			return values.size();
		};
	};
}