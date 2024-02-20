#include "Square.hpp"
#include "tanh_sinh_helper.hpp"
#include "../Constants.hpp"
#include <numeric>
#include <filesystem>

namespace Hubbard::DensityOfStates {
	const std::string DATA_FILE_NAME{ "../../data/dos_square.bin" };

	std::vector<abscissa_t> Square::upper_border_to_abscissa;
	dos_precision Square::b_minus_a_halved;

	global_floating_type Square::computeValue(const global_floating_type& gamma) {
		return (LONG_1_PI * LONG_1_PI) * boost::math::ellint_1(sqrt(1 - 0.25 * gamma * gamma));
	};

#ifdef _BOOST_PRECISION
#pragma omp declare reduction(+:dos_precision:omp_out+=omp_in)
#endif
	void Square::computeValues()
	{
		clearAll();
		if (std::filesystem::exists(DATA_FILE_NAME)) {
			if (loadFromBinaryFile(DATA_FILE_NAME)) {
				return;
			}
			else {
				// We end up here if something went wrong while loading
				// In this case we want a fresh computation
				clearAll();
			}
		}

		step = std::ldexp(1, -1);
		auto compute_DOS = [](abscissa_t gamma, abscissa_t one_minus_gamma) -> dos_precision {
			gamma *= 0.5;
			one_minus_gamma *= 0.5;
			return (LONG_1_PI * LONG_1_PI)
				* boost::math::ellint_1((gamma < 0.25 ? sqrt(1 - gamma * gamma) : sqrt(one_minus_gamma * (1 + gamma)))).convert_to<dos_precision>();
			};

		tanh_sinh_helper<abscissa_t, dos_precision> tsh{ 0, 2 };
		tanh_sinh_helper<abscissa_t, dos_precision>::SaveTo buffer_vectors{ &abscissa, &upper_border_to_abscissa, &weights, &values };
		dos_precision old_integral{ tsh.initial_filling<SQUARE_QUAD_CUT_OFF>(compute_DOS, buffer_vectors) };

		dos_precision new_integral{};
		dos_precision error{ 100.0 };
		while (error > (boost::math::pow<SQUARE_QUAD_CUT_OFF>(10))) {
			tsh.increase_level(buffer_vectors);

			new_integral = 0;
#pragma omp parallel for reduction(+:new_integral) schedule(dynamic)
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

			error = abs(new_integral - old_integral);
			std::cout << error << std::endl;
			old_integral = new_integral;
		}

		std::cout << "Square: Exit after " << tsh.level() << " levels with error = " << abs(0.5 - new_integral) << std::endl;
		std::cout << "Total amount of values = " << values.size() << std::endl;

		computed = true;
		writeToBinaryFile(DATA_FILE_NAME);
	}
}