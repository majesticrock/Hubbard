#define _USE_MATH_DEFINES

#include "Square.hpp"
#include "../Constants.hpp"
#include <boost/math/special_functions/ellint_1.hpp>

namespace Hubbard::DensityOfStates {
	void Square::computeValues()
	{
		step = 2. / Constants::K_DISCRETIZATION;
		values.resize(Constants::K_DISCRETIZATION);
		for (int g = -Constants::K_DISCRETIZATION; g < 0; g++)
		{
			const double gamma = (0.5 + g) * step;
			values[g + Constants::K_DISCRETIZATION] = (M_1_PI * M_1_PI)
				* boost::math::ellint_1(sqrt(1. - (gamma * gamma / 4.)));
		}
		computed = true;
	}
}