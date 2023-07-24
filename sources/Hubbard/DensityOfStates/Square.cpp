#define _USE_MATH_DEFINES

#include "Square.hpp"
#include "../Constants.hpp"
#include <boost/math/special_functions/ellint_1.hpp>

namespace Hubbard::DensityOfStates {
	void Square::computeValues()
	{
		step = 2. / Constants::BASIS_SIZE;
		values.reserve(2 * Constants::BASIS_SIZE);
		values.resize(Constants::BASIS_SIZE);

		for (int g = 0; g < Constants::BASIS_SIZE; ++g)
		{
			const double gamma = (0.5 + g) * step;
			values[g] = (M_1_PI * M_1_PI)
				* boost::math::ellint_1(sqrt(1. - (gamma * gamma / 4.)));
		}
		symmetrizeVector<false>(values);
		computed = true;
	}
}