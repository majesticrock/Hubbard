#define _USE_MATH_DEFINES

#include "Square.hpp"
#include <boost/math/special_functions/ellint_1.hpp>
#include <cmath>
#include "../Constants.hpp"

namespace Hubbard::DensityOfStates {
	std::vector<double> Square::values;
	bool Square::computed = false;
	double Square::step = 0;

	void Square::computeValues()
	{
		Square::step = 2. / Constants::K_DISCRETIZATION;
		Square::values.resize(Constants::K_DISCRETIZATION);
		for (int g = -Constants::K_DISCRETIZATION; g < 0; g++)
		{
			const double gamma = (0.5 + g) * Square::step;
			Square::values[g + Constants::K_DISCRETIZATION] = (M_1_PI * M_1_PI)
				* boost::math::ellint_1(sqrt(1. - (gamma * gamma / 4.)));
		}
		Square::computed = true;
	}
}