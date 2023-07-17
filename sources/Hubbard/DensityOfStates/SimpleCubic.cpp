#include "SimpleCubic.hpp"
#include "../Constants.hpp"
#include <boost/math/special_functions/ellint_1.hpp>

namespace Hubbard::DensityOfStates {
	void SimpleCubic::computeValues()
	{
		step = 3. / Constants::BASIS_SIZE;
		values.resize(Constants::BASIS_SIZE);

	}
}