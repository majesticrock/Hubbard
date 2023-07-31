#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

namespace Hubbard::DensityOfStates {
	constexpr long double LONG_PI = 3.141592653589793238462643383L;
	constexpr long double LONG_1_PI = 1 / LONG_PI;
	constexpr long double LOG_4 = 2 * 0.693147180559945309417L; // ln(4) = 2 ln(2)
	constexpr long double R_AT_2 = (LONG_PI - 4.L) / 8;
	constexpr long double R_AT_2_1 = (5 * LONG_PI / 64 - 0.25L);
	constexpr long double CUT_OFF = 1e-12;

	// Computes sqrt(1 - x^2)
	template <class RealType>
	inline RealType sqrt_1_minus_x_squared(RealType x) {
		return std::sqrt(1 - x * x);
	}

	inline long double R(long double x) {
		x *= 0.5;
		// Special, analytically known cases:
		if (std::abs(x) < CUT_OFF) {
			return LOG_4;
		}
		return boost::math::ellint_1(sqrt_1_minus_x_squared(x)) + 0.5 * std::log(x * x);
	}

	inline long double derivative_R(long double x) {
		// Special, analytically known cases:
		if (std::abs(x) < CUT_OFF) {
			return 0.0L;
		}
		if (std::abs(x + 2) < CUT_OFF) {
			return R_AT_2 + R_AT_2_1 * (x + 2);
		}
		if (std::abs(x - 2) < CUT_OFF) {
			return -R_AT_2 + R_AT_2_1 * (x - 2);
		}

		const long double ALPHA = sqrt_1_minus_x_squared(0.5 * x);
		return (x * x * (1 - boost::math::ellint_1(ALPHA)) + 4 * (boost::math::ellint_2(ALPHA) - 1)) / (x * (x + 2) * (x - 2));
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
	inline double computeNorm() {
		return typename DOS::DOSIntegrator<double>().integrate_by_value([](double) -> double { return 1.0; });
	}

	struct BaseDOS {
		// Contains the values for gamma in (-2, 0) as the the positive part is symmetric.
		// Open boundaries, as the dos at 0 contains a singularity for d=2
		static std::vector<double> values;
		static double step;
		static bool computed;
		static double integrateValues();
		static void printValues();

		virtual void computeValues() = 0;
		virtual ~BaseDOS() = default;
	};
}