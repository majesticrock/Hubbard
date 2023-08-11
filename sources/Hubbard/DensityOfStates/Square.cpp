#include "Square.hpp"
#include "tanh_sinh_helper.hpp"
#include "../Constants.hpp"
#include <numeric>

namespace Hubbard::DensityOfStates {
	std::vector<abscissa_t> Square::abscissa;
	std::vector<abscissa_t> Square::upper_border_to_abscissa;
	std::vector<double> Square::weights;

	double Square::LOWER_BORDER;
	double Square::b_minus_a_halved;

	void Square::computeValues()
	{
		step = std::ldexp(1, -1);
		auto compute_DOS = [](abscissa_t gamma, abscissa_t one_minus_gamma) -> double {
			gamma *= 0.5;
			one_minus_gamma *= 0.5;
			return (LONG_1_PI * LONG_1_PI)
				* boost::math::ellint_1((gamma < 0.25 ? sqrt(1 - gamma * gamma) : sqrt(one_minus_gamma * (1 + gamma)))).convert_to<double>();
		};

		tanh_sinh_helper<abscissa_t> tsh{ 0, 2 };
		tanh_sinh_helper<abscissa_t>::SaveTo buffer_vectors{ &abscissa, & upper_border_to_abscissa, & weights, & values};
		double old_integral{ tsh.initial_filling(compute_DOS, buffer_vectors) };

		double new_integral{};
		double error{ 100.0 };
		while (error > 1e-12) {
			tsh.increase_level(buffer_vectors);

			new_integral = 0;
			for (int k = 0; k < values.size(); ++k)
			{
				if (k % 2 != 0) {
					tsh.compute_step(compute_DOS, k, buffer_vectors);
				}
				else {
					weights[k] *= 0.5;
				}
				new_integral += values[k] * weights[k];
			}
			new_integral *= tsh.half_distance();

			error = std::abs(new_integral - old_integral);
			old_integral = new_integral;
		}

		std::cout << "Exit after " << tsh.level() << " levels with error = " << std::abs(0.5 - new_integral) << std::endl;
		std::cout << "Total amount of values = " << values.size() << std::endl;

		computed = true;
	}
}