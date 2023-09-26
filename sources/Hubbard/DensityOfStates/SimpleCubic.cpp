#include "SimpleCubic.hpp"
#include "../Constants.hpp"
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "tanh_sinh_helper.hpp"

namespace Hubbard::DensityOfStates {
	using std::log;
	using std::asin;

	std::vector<abscissa_t> SimpleCubic::upper_border_to_abscissa;
	dos_precision SimpleCubic::b_minus_a_halved;

	using _internal_precision = dos_precision;

	inline _internal_precision R(abscissa_t x) {
		x *= 0.5;
		// Special, analytically known cases:
		if (abs(x) < CUT_OFF) {
			return LOG_4;
		}
		return static_cast<_internal_precision>(boost::math::ellint_1(sqrt_1_minus_x_squared(x)) + 0.5 * log(x * x));
	}

	inline _internal_precision derivative_R(abscissa_t x) {
		// Special, analytically known cases:
		if (abs(x) < CUT_OFF) {
			return 0.0;
		}
		if (abs(x + 2) < CUT_OFF) {
			return R_AT_2 + R_AT_2_1 * static_cast<_internal_precision>(x + 2);
		}
		if (abs(x - 2) < CUT_OFF) {
			return -R_AT_2 + R_AT_2_1 * static_cast<_internal_precision>(x - 2);
		}

		const abscissa_t ALPHA = sqrt_1_minus_x_squared(0.5 * x);
		return static_cast<_internal_precision>((x * x * (1 - boost::math::ellint_1(ALPHA)) + 4 * (boost::math::ellint_2(ALPHA) - 1)) / (x * (x + 2) * (x - 2)));
	}

	inline _internal_precision I_1(abscissa_t gamma) {
		const abscissa_t lower_bound{ std::max(abscissa_t{ -1 }, abscissa_t{ -2 } - gamma) };
		const abscissa_t upper_bound{ std::min(abscissa_t{ 1 }, abscissa_t{ 2 } - gamma) };
		_internal_precision ret = static_cast<_internal_precision>(asin(upper_bound) * R(upper_bound + gamma) - asin(lower_bound) * R(lower_bound + gamma));

		auto integrand = [gamma](abscissa_t phi) {
			return asin(phi) * derivative_R(phi + gamma);
			};

		ret -= static_cast<_internal_precision>(boost::math::quadrature::gauss_kronrod<abscissa_t, 30>::integrate(integrand, lower_bound, upper_bound, 10, 1e-12));
		return ret;

		//boost::math::quadrature::tanh_sinh<abscissa_t> integrator;
		//ret -= static_cast<_internal_precision>(integrator.integrate(integrand, lower_bound, upper_bound));
		return ret;
	}

	inline _internal_precision I_2(abscissa_t gamma) {
		// For some magical reason this integral is constant for gamma in [-1, 1]
		// I_2 = -pi * ln(4)
		if (gamma >= -1 && gamma <= 1) {
			return (-LONG_PI * LOG_4);
		}
		const abscissa_t lower_bound{ std::max(abscissa_t{ -1 }, abscissa_t{ -2 } - gamma) };
		const abscissa_t upper_bound{ std::min(abscissa_t{ 1 }, abscissa_t{ 2 } - gamma) };

		//std::cout << "Bounds: " << gamma << "   -   " << lower_bound << "  " << upper_bound << std::endl;

		auto integrand = [gamma](abscissa_t phi) {
			return log(0.25 * (gamma + phi) * (gamma + phi)) / sqrt_1_minus_x_squared(phi);
			};
		boost::math::quadrature::tanh_sinh<abscissa_t> integrator;
		return 0.5 * static_cast<_internal_precision>(integrator.integrate(integrand, lower_bound, upper_bound));
	}

	template <class T>
	inline auto append_vectors(std::vector<T>& a, const std::vector<T>& b) {
		return a.insert(a.end(), b.begin(), b.end());
	}

#ifdef _BOOST_PRECISION
#pragma omp declare reduction(+:dos_precision:omp_out+=omp_in)
#endif
	void SimpleCubic::computeValues()
	{
		step = std::ldexp(1, -1);
		auto compute_DOS = [](abscissa_t gamma, abscissa_t one_minus_gamma) -> dos_precision {
			return static_cast<dos_precision>(boost::math::pow<3>(LONG_1_PI) * (I_1(gamma) - I_2(gamma)));
			};

		constexpr double borders[2][2] = { {0, 1},{1, 3} };

		for (size_t i = 0U; i < 2U; ++i)
		{
			decltype(abscissa) buf_abscissa;
			decltype(upper_border_to_abscissa) buf_upper_border_to_abscissa;
			decltype(weights) buf_weights;
			decltype(values) buf_values;

			tanh_sinh_helper<abscissa_t, dos_precision> tsh{ borders[i][0], borders[i][1] };
			tanh_sinh_helper<abscissa_t, dos_precision>::SaveTo buffer_vectors{ &buf_abscissa, &buf_upper_border_to_abscissa, &buf_weights, &buf_values };
			dos_precision old_integral{ tsh.initial_filling<SC_QUAD_CUT_OFF>(compute_DOS, buffer_vectors) };

			dos_precision new_integral{};
			dos_precision error{ 100.0 };
			while (error > (boost::math::pow<SC_QUAD_CUT_OFF>(10))) {
				tsh.increase_level(buffer_vectors);

				new_integral = 0;
#pragma omp parallel for reduction(+:new_integral) schedule(dynamic)
				for (int k = 0; k < buf_values.size(); ++k)
				{
					if (k % 2 != 0) {
						tsh.compute_step(compute_DOS, k, buffer_vectors);
					}
					else {
						buf_weights[k] *= 0.5;
					}
					new_integral += buf_values[k] * buf_weights[k];
				}
				new_integral *= tsh.half_distance();

				error = abs(new_integral - old_integral);
				std::cout << error << std::endl;
				old_integral = new_integral;
			}

			std::cout << "Exit after " << tsh.level() << " levels with error = " << error << std::endl;
			std::cout << "Total amount of values = " << buf_values.size() << std::endl;
			if (i == 0) {
				abscissa = std::move(buf_abscissa);
				upper_border_to_abscissa = std::move(buf_upper_border_to_abscissa);
				weights = std::move(buf_weights);
				values = std::move(buf_values);
			}
			else {
				append_vectors(abscissa, buf_abscissa);
				append_vectors(upper_border_to_abscissa, buf_upper_border_to_abscissa);
				append_vectors(weights, buf_weights);
				append_vectors(values, buf_values);
			}
		}

		computed = true;
	}
}