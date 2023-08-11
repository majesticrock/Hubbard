#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

namespace Hubbard::DensityOfStates {
	typedef long double dos_precision;
	constexpr long double LONG_PI = 3.141592653589793238462643383L;
	constexpr long double LONG_PI_2 = LONG_PI / 2;
	constexpr long double LONG_1_PI = 1 / LONG_PI;
	constexpr long double LOG_4 = 2 * 0.693147180559945309417L; // ln(4) = 2 ln(2)
	constexpr long double R_AT_2 = (LONG_PI - 4.L) / 8;
	constexpr long double R_AT_2_1 = (5 * LONG_PI / 64 - 0.25L);
	constexpr long double CUT_OFF = 1e-12;

	// Computes sqrt(1 - x^2)
	template <class RealType>
	inline RealType sqrt_1_minus_x_squared(RealType x) {
		return std::sqrt((1 - x) * (1 + x));
	}
	// Use this overload for better results, if |x| is close to 1
	// However, you need to provide an accurate reprensantion of 1+/-x yourself.
	template <class RealType>
	inline RealType sqrt_1_minus_x_squared(RealType x, RealType one_minus_x) {
		return (x < 0.25 ? std::sqrt(1 - x * x) : std::sqrt(one_minus_x * (1 + x)));
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
		return typename DOS::DOSIntegrator<dos_precision>().integrate_by_value([](dos_precision) -> dos_precision { return 1.L; });
	}

	struct BaseDOS {
		// Contains the values for gamma in (-2, 0) as the the positive part is symmetric.
		// Open boundaries, as the dos at 0 contains a singularity for d=2
		static std::vector<dos_precision> values;
		static dos_precision step;
		static bool computed;
		static dos_precision integrateValues();
		static void printValues();

		virtual void computeValues() = 0;
		virtual ~BaseDOS() = default;
	};
}