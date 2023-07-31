#include "SimpleCubic.hpp"
#include "../Constants.hpp"
#include <omp.h>
#include <cmath>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace Hubbard::DensityOfStates {
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
		constexpr size_t num_positions = 30U;
		step = 3. / Constants::BASIS_SIZE;
		const int splits_s[] = { 3, 6, 12, 24, 48, 96, 192, 384 };
		for (const int splits : splits_s) {
			values.clear();
			values.resize(splits * num_positions / 2);
			/*
			*  Algorithm, remember, we have t = 0.5
			*/
			auto getIndex = [](double k) {
				for (size_t i = 0U; i < num_positions; ++i)
				{
					if (std::abs(boost::math::quadrature::gauss<double, num_positions>::abscissa()[i] - k) < 1e-12
						|| std::abs(boost::math::quadrature::gauss<double, num_positions>::abscissa()[i] + k) < 1e-12) {
						return i;
					}
				}
				throw std::runtime_error("Did not find k within the abscissa!");
			};

			int offset = 0;
			std::pair<double, double> lims{-3, -1};

			auto procedure = [&](long double gamma) -> double {
				double k = (gamma - 0.5 * (lims.first + lims.second)) * (2. / (lims.second - lims.first));
				const int pos = offset + getIndex(k);
				values[pos] = boost::math::pow<3>(M_1_PI) * (I_1(gamma) - I_2(gamma));
				return values[pos];
			};

			double norm = 0;
			for (int i = 0; i < splits; ++i)
			{
				lims = { -3. + i * 6. / splits, -3. + (i + 1) * 6. / splits };
				norm += boost::math::quadrature::gauss<double, num_positions>::integrate(procedure, lims.first, lims.second);
				offset += num_positions / 2;
			}

			//symmetrizeVector<true>(values);
			std::cout << splits << ":    1 - NORM = " << std::scientific << 1. - norm << std::endl;
		}
		computed = true;
	}
}