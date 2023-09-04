#include "SimpleCubic.hpp"
#include "../Constants.hpp"
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace Hubbard::DensityOfStates {
	using std::log;
	using std::asin;

	std::vector<std::pair<dos_precision, dos_precision>> SimpleCubic::split_limits;
	std::array<dos_precision, SimpleCubic::num_positions> SimpleCubic::abscissa;
	std::array<dos_precision, SimpleCubic::num_positions> SimpleCubic::weights;
	int SimpleCubic::n_splits;
	dos_precision SimpleCubic::b_minus_a_halved;

	using _internal_precision = dos_precision;

	inline _internal_precision R(_internal_precision x) {
		x *= 0.5;
		// Special, analytically known cases:
		if (abs(x) < CUT_OFF) {
			return LOG_4;
		}
		return boost::math::ellint_1(sqrt_1_minus_x_squared(x)) + 0.5 * log(x * x);
	}

	inline _internal_precision derivative_R(_internal_precision x) {
		// Special, analytically known cases:
		if (abs(x) < CUT_OFF) {
			return 0.0L;
		}
		if (abs(x + 2) < CUT_OFF) {
			return R_AT_2 + R_AT_2_1 * (x + 2);
		}
		if (abs(x - 2) < CUT_OFF) {
			return -R_AT_2 + R_AT_2_1 * (x - 2);
		}

		const _internal_precision ALPHA = sqrt_1_minus_x_squared(0.5 * x);
		return (x * x * (1 - boost::math::ellint_1(ALPHA)) + 4 * (boost::math::ellint_2(ALPHA) - 1)) / (x * (x + 2) * (x - 2));
	}

	inline _internal_precision I_1(_internal_precision gamma) {
		const _internal_precision lower_bound = std::max(_internal_precision{ -1 }, _internal_precision{ -2 } - gamma);
		const _internal_precision upper_bound = std::min(_internal_precision{ 1 }, _internal_precision{ 2 } - gamma);
		_internal_precision ret = asin(upper_bound) * R(upper_bound + gamma) - asin(lower_bound) * R(lower_bound + gamma);

		auto integrand = [gamma](_internal_precision phi) {
			return asin(phi) * derivative_R(phi + gamma);
		};

		ret -= boost::math::quadrature::gauss_kronrod<_internal_precision, 30>::integrate(integrand, lower_bound, upper_bound, 10, 1e-12);
		return ret;
	}

	inline _internal_precision I_2(_internal_precision gamma) {
		// For some magical reason this integral is constant for gamma in [-1, 1]
		// I_2 = -pi * ln(4)
		if (gamma >= -1 && gamma <= 1) {
			return (-LONG_PI * LOG_4);
		}
		const _internal_precision lower_bound = std::max(_internal_precision{ -1 }, _internal_precision{ -2 } - gamma);
		const _internal_precision upper_bound = std::min(_internal_precision{ 1 }, _internal_precision{ 2 } - gamma);

		auto integrand = [gamma](_internal_precision phi) {
			return log(0.25 * (gamma + phi) * (gamma + phi)) / sqrt_1_minus_x_squared(phi);
		};
		boost::math::quadrature::tanh_sinh<_internal_precision> integrator;
		return 0.5 * integrator.integrate(integrand, lower_bound, upper_bound);
	}

	void SimpleCubic::computeValues()
	{
		n_splits = 3 * Constants::K_DISCRETIZATION;
		b_minus_a_halved = 1. / Constants::K_DISCRETIZATION;

		step = 3. / Constants::BASIS_SIZE;
		values.clear();
		values.resize(n_splits * num_positions);
		split_limits.resize(n_splits);

		std::copy(boost::math::quadrature::gauss<dos_precision, num_positions>::abscissa().begin(),
			boost::math::quadrature::gauss<dos_precision, num_positions>::abscissa().end(), abscissa.begin());
		for (auto& a : abscissa)
		{
			a *= -1;
		}
		std::copy(boost::math::quadrature::gauss<dos_precision, num_positions>::abscissa().rbegin(),
			boost::math::quadrature::gauss<dos_precision, num_positions>::abscissa().rend(), abscissa.begin() + num_positions / 2);
		for (auto& a : abscissa)
		{
			a *= b_minus_a_halved;
		}
		std::cout << std::endl;
		std::copy(boost::math::quadrature::gauss<dos_precision, num_positions>::weights().begin(),
			boost::math::quadrature::gauss<dos_precision, num_positions>::weights().end(), weights.begin());
		std::copy(boost::math::quadrature::gauss<dos_precision, num_positions>::weights().rbegin(),
			boost::math::quadrature::gauss<dos_precision, num_positions>::weights().rend(), weights.begin() + num_positions / 2);
		for (auto& w : weights)
		{
			w *= b_minus_a_halved;
		}

		/*
		*  Algorithm, remember, we have t = 0.5
		*/
		for (int i = 0; i < n_splits; ++i)
		{
			split_limits[i] = { -3. + i * 6. / n_splits, -3. + (i + 1) * 6. / n_splits };
		}
		for (int i = 0; i < n_splits; ++i)
		{
#pragma omp parallel for
			for (size_t j = 0U; j < num_positions; ++j)
			{
				const int pos = i * num_positions + j;
				const dos_precision gamma = functionQuadratureOffset(i) + abscissa[j];

				values[pos] = boost::math::pow<3>(LONG_1_PI) * (I_1(gamma) - I_2(gamma));
			}
		}

		computed = true;
	}
}