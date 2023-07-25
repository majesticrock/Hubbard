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
		step = 3. / Constants::BASIS_SIZE;
		const int splits_s[] = { 3, 6, 12, 24, 48, 96, 192, 384 };
		for (const int splits : splits_s) {
			values.clear();
			values.resize(splits * num_positions / 2);
			/*
			*  Algorithm, remember, we have t = 0.5
			*/
			boost::math::quadrature::tanh_sinh<long double> integrator;
			constexpr double singularity_offset = 1e-7;

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

			auto procedure = [&](double gamma) -> double {
				const double lower_bound = std::max(-1., -2. - gamma);
				const double upper_bound = std::min(1., 2. - gamma);
				auto boost_integrand = [gamma](long double phi) -> long double {
					return boost::math::ellint_1(sqrt_1_minus_x_squared((gamma + phi) / 2.L)) / sqrt_1_minus_x_squared(phi);
				};

				double k = (gamma - 0.5 * (lims.first + lims.second)) * (2. / (lims.second - lims.first));
				const int pos = offset + getIndex(k);
				if (gamma <= -3 || gamma >= 3) {
					values[pos] = 0.;
					return 0.;
				}
				const long double tol = std::sqrt(std::numeric_limits<double>::epsilon());
				long double error = 0;
				long double L1 = 0;
				size_t levels = 0;
				if (-gamma - singularity_offset > lower_bound && -gamma + singularity_offset < upper_bound) {
					values[pos] = boost::math::pow<3>(M_1_PI)
						* integrator.integrate(boost_integrand, lower_bound, -gamma - singularity_offset, tol, &error, &L1, &levels);
					values[pos] += boost::math::pow<3>(M_1_PI)
						* integrator.integrate(boost_integrand, -gamma + singularity_offset, upper_bound, tol, &error, &L1, &levels);
				}
				else {
					values[pos] = boost::math::pow<3>(M_1_PI)
						* integrator.integrate(boost_integrand, lower_bound, upper_bound, tol, &error, &L1, &levels);
				}
				std::cout << error << "   " << L1 << "   " << levels << std::endl;
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