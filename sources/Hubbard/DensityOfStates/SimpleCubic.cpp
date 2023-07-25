#define _USE_MATH_DEFINES
#include "SimpleCubic.hpp"
#include "../Constants.hpp"
#include <omp.h>
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace Hubbard::DensityOfStates {
	constexpr long double LONG_PI = 3.141592653589793238462643383L;
	constexpr long double LOG_4 = 2 * 0.693147180559945309417L; // ln(4) = 2 ln(2)
	constexpr long double R_AT_2 = (LONG_PI - 4.L) / 8;
	constexpr long double R_AT_2_1 = (5 * LONG_PI / 64 - 0.25L);
	constexpr long double CUT_OFF = 1e-12;

	using boost::math::ellint_1;
	using boost::math::ellint_2;

	inline long double R(long double phi, long double gamma) {
		// Special, analytically known cases:
		const long double SUM = 0.5 * (phi + gamma);
		if (std::abs(SUM) < CUT_OFF) {
			return LOG_4;
		}
		return ellint_1(sqrt_1_minus_x_squared(SUM)) + 0.5 * std::log(SUM * SUM);
	}

	inline long double derivative_R(long double phi, long double gamma) {
		// Special, analytically known cases:
		const long double SUM = phi + gamma;
		if (std::abs(SUM) < CUT_OFF) {
			return 0.0L;
		}
		if (std::abs(SUM + 2) < CUT_OFF) {
			return R_AT_2 * R_AT_2_1 * (SUM + 2);
		}
		if (std::abs(SUM - 2) < CUT_OFF) {
			return -R_AT_2 + R_AT_2_1 * (SUM - 2);
		}

		const long double ALPHA = sqrt_1_minus_x_squared(0.5 * SUM);
		return (SUM * SUM * (1 - ellint_1(ALPHA)) + 4 * (ellint_2(ALPHA) - 1)) / (SUM * (SUM + 2) * (SUM - 2));
	}

	inline long double I_1(long double gamma) {
		const long double lower_bound = std::max(-1.L, -2.L - gamma);
		const long double upper_bound = std::min(1.L, 2.L - gamma);
		long double ret = std::asin(upper_bound) * R(upper_bound, gamma) - std::asin(lower_bound) * R(lower_bound, gamma);

		auto integrand = [gamma](long double phi) {
			return std::asin(phi) * derivative_R(phi, gamma);
		};
		ret -= boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, lower_bound, upper_bound, 10, 1e-12);
		return ret;
	}

	inline long double I_2(long double gamma) {
		// For some magical reason this integral is constant for gamma in [-1, 1]
		// I_2 = -pi * ln(4)
		if (gamma >= -1 && gamma <= 1) {
			return (- LONG_PI * LOG_4);
		}
		const long double lower_bound = std::max(-1.L, -2.L - gamma);
		const long double upper_bound = std::min(1.L, 2.L - gamma);

		auto integrand = [gamma](long double phi) {
			return std::log(0.25 * (gamma + phi) * (gamma + phi)) / sqrt_1_minus_x_squared(phi);
		};
		boost::math::quadrature::tanh_sinh<double> integrator;
		return 0.5 * integrator.integrate(integrand, lower_bound, upper_bound);
	}

	void SimpleCubic::computeValues()
	{
		step = 3. / Constants::BASIS_SIZE;
		values.reserve(2 * Constants::BASIS_SIZE - 1);
		values.resize(Constants::BASIS_SIZE + 1);
		values.back() = 0.0;

#pragma omp parallel for
		for (int g = 0; g < Constants::BASIS_SIZE; ++g)
		{
			const double gamma = g * step;
			values[g] = boost::math::pow<3>(M_1_PI) * (I_1(gamma) - I_2(gamma));
		}
		symmetrizeVector<true>(values);
		computed = true;
	}
}