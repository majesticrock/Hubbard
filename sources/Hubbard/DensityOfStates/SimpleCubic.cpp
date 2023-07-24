#define _USE_MATH_DEFINES
#include "SimpleCubic.hpp"
#include "../Constants.hpp"
#include <omp.h>
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp>

namespace Hubbard::DensityOfStates {
	// Computes sqrt(1 - x^2)
	template <class RealType>
	inline RealType sqrt_1_minus_x_squared(RealType x) {
		return sqrt(1 - x * x);
	}

	void SimpleCubic::computeValues()
	{
		constexpr size_t num_positions = 30U;
		const auto& abcissa = boost::math::quadrature::gauss<double, num_positions>::abscissa();

		step = 3. / Constants::BASIS_SIZE;
		values.reserve(2 * Constants::BASIS_SIZE - 1);
		values.resize(Constants::BASIS_SIZE + 1);
		values.back() = 0.0;
		/*
		*  Algorithm, remember, we have t = 0.5
		*/
		boost::math::quadrature::tanh_sinh<double> integrator;
		constexpr double singularity_offset = 1e-7;

#pragma omp parallel for
		for (int g = 0; g < Constants::BASIS_SIZE; ++g)
		{
			const double gamma = (g)*step;
			constexpr double lower_bound = -1;// std::max(-1., -2. - gamma);
			const double upper_bound = std::min(1., 2. - gamma);

			auto boost_integrand = [gamma](long double phi) -> double {
				return boost::math::ellint_1(sqrt_1_minus_x_squared((gamma + phi) / 2.L)) / sqrt_1_minus_x_squared(phi);
			};

			if (-gamma - singularity_offset > lower_bound) {
				values[g] = boost::math::pow<3>(M_1_PI)
					* integrator.integrate(boost_integrand, lower_bound, -gamma - singularity_offset);
				values[g] += boost::math::pow<3>(M_1_PI)
					* integrator.integrate(boost_integrand, -gamma + singularity_offset, upper_bound);
			}
			else {
				values[g] = boost::math::pow<3>(M_1_PI)
					* integrator.integrate(boost_integrand, lower_bound, upper_bound);
			}
		}
		symmetrizeVector<true>(values);
		computed = true;
	}
}