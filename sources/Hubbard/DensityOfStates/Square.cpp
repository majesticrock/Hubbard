#include "Square.hpp"
#include "../Constants.hpp"
#include <numeric>

namespace Hubbard::DensityOfStates {
	std::vector<double> Square::regular_values;
	std::vector<double> Square::singular_values_linear;
	std::vector<double> Square::singular_values_quadratic;

	std::vector<double> Square::singular_weights;
	double Square::LOWER_BORDER;

	template <class RealType>
	// Returns x*ln((x/2)^2)
	inline RealType x_log_x_2_squared(RealType x) {
		if (std::abs(x) < CUT_OFF) return 0;
		return x * (std::log(x * x) - LOG_4);
	}

	void Square::computeValues()
	{
		step = 2. / Constants::BASIS_SIZE;
		LOWER_BORDER = -2 + 0.5 * step;

		values.reserve(2 * Constants::BASIS_SIZE);
		values.resize(Constants::BASIS_SIZE);

		for (int g = 0; g < Constants::BASIS_SIZE; ++g)
		{
			const long double gamma = (0.5 + g) * step;
			values[g] = (LONG_1_PI * LONG_1_PI)
				* boost::math::ellint_1(sqrt_1_minus_x_squared(0.5 * gamma));
		}
		symmetrizeVector<false>(values);

		const size_t VECTOR_SIZE = 2 * Constants::BASIS_SIZE;
		regular_values.resize(VECTOR_SIZE);
		singular_values_linear.resize(VECTOR_SIZE);
		singular_values_quadratic.resize(VECTOR_SIZE);

		for (int g = 0; g < VECTOR_SIZE; ++g)
		{
			const long double gamma = LOWER_BORDER + g * step;
			regular_values[g] = (LONG_1_PI * LONG_1_PI) * R(gamma);
			singular_values_linear[g] = (LONG_1_PI * LONG_1_PI) * (x_log_x_2_squared(gamma) - 2 * gamma);
			singular_values_quadratic[g] = (LONG_1_PI * LONG_1_PI) * 0.5 * gamma * (x_log_x_2_squared(gamma) - gamma);
		}

		auto weight = [](double gamma){
			if(gamma > -1e-12 && gamma < 1e-12) return 0.0;
			return gamma * (std::log(gamma * gamma * 0.25) - 2);
		};

		singular_weights.resize(VECTOR_SIZE);

		for (int g = 1; g < VECTOR_SIZE; g++)
		{
			long double gamma = LOWER_BORDER + (g - 0.5) * step;
			singular_weights[g] = LONG_1_PI * LONG_1_PI * (weight(gamma + step) - weight(gamma));
		}
		
		computed = true;
	}
}