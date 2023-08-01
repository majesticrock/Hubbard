#include "SimpleCubic.hpp"
#include "../Constants.hpp"
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace Hubbard::DensityOfStates {
	std::vector<std::pair<double, double>> SimpleCubic::split_limits;
	std::array<double, SimpleCubic::num_positions> SimpleCubic::abscissa;
	std::array<double, SimpleCubic::num_positions> SimpleCubic::weights;
	int SimpleCubic::n_splits;
	double SimpleCubic::b_minus_a_halved;

	inline long double I_1(long double gamma) {
		const long double lower_bound = std::max(-1.L, -2.L - gamma);
		const long double upper_bound = std::min(1.L, 2.L - gamma);
		long double ret = std::asin(upper_bound) * R(upper_bound + gamma) - std::asin(lower_bound) * R(lower_bound + gamma);

		auto integrand = [gamma](long double phi) {
			return std::asin(phi) * derivative_R(phi + gamma);
		};
		ret -= boost::math::quadrature::gauss_kronrod<double, 15>::integrate(integrand, lower_bound, upper_bound, 10, 1e-12);
		return ret;
	}

	inline long double I_2(long double gamma) {
		// For some magical reason this integral is constant for gamma in [-1, 1]
		// I_2 = -pi * ln(4)
		if (gamma >= -1 && gamma <= 1) {
			return (-LONG_PI * LOG_4);
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
		n_splits = 3 * Constants::K_DISCRETIZATION;
		b_minus_a_halved = 1. / Constants::K_DISCRETIZATION;

		step = 3. / Constants::BASIS_SIZE;
		values.clear();
		values.resize(n_splits * num_positions);
		split_limits.resize(n_splits);

		std::copy(boost::math::quadrature::gauss<double, num_positions>::abscissa().begin(),
			boost::math::quadrature::gauss<double, num_positions>::abscissa().end(), abscissa.begin());
		for (auto& a : abscissa)
		{
			a *= -1;
		}
		std::copy(boost::math::quadrature::gauss<double, num_positions>::abscissa().rbegin(),
			boost::math::quadrature::gauss<double, num_positions>::abscissa().rend(), abscissa.begin() + num_positions / 2);
		for (auto& a : abscissa)
		{
			a *= b_minus_a_halved;
		}
		std::cout << std::endl;
		std::copy(boost::math::quadrature::gauss<double, num_positions>::weights().begin(),
			boost::math::quadrature::gauss<double, num_positions>::weights().end(), weights.begin());
		std::copy(boost::math::quadrature::gauss<double, num_positions>::weights().rbegin(),
			boost::math::quadrature::gauss<double, num_positions>::weights().rend(), weights.begin() + num_positions / 2);
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
			for (size_t j = 0U; j < num_positions; ++j)
			{
				const int pos = i * num_positions + j;
				const double gamma = functionQuadratureOffset(i) + abscissa[j];

				values[pos] = boost::math::pow<3>(LONG_1_PI) * (I_1(gamma) - I_2(gamma));
			}
		}

		computed = true;
	}
}