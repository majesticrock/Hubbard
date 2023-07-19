#define _USE_MATH_DEFINES
#include "SimpleCubic.hpp"
#include "../Constants.hpp"
#include "../../Utility/MidpointRule.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/pow.hpp>

namespace Hubbard::DensityOfStates {
	void SimpleCubic::computeValues()
	{
		constexpr size_t MID_POINT_STEP_NUM = 10000;
		/*
		*  Algorithm, remember, we have t = 0.5
		*/
		step = 3. / Constants::BASIS_SIZE;
		values.resize(Constants::BASIS_SIZE);

		auto integrand = [](const double phi, const double gamma) -> double {
			const double pass_to_rho_2D = gamma + phi;
			return boost::math::ellint_1(sqrt(1. - (pass_to_rho_2D * pass_to_rho_2D / 4.))) / sqrt(1 - phi * phi);
		};

		for (int g = 0; g < Constants::BASIS_SIZE; ++g)
		{
			const double gamma = (0.5 + g) * step;
			const double lower_bound = std::max(-1., -2. - gamma);
			const double upper_bound = std::min(1., 2. - gamma);
			// We extractet 1/pi^2 from the 2D DOS for computational benefit
			values[g] = boost::math::pow<3>(M_1_PI) * Utility::NumericalSolver::Integration::midpoint_rule(integrand, lower_bound, upper_bound, MID_POINT_STEP_NUM, gamma);
		}
		computed = true;

		// Factor of 2 because we only computed half the DOS due to symmetry
		std::cout << "Norm: " << 2 * std::accumulate(values.begin(), values.end(), double{}) * step << std::endl;
	}
}