#define _USE_MATH_DEFINES
#include "SimpleCubic.hpp"
#include "../Constants.hpp"
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace Hubbard::DensityOfStates {
	template <class UnaryFunction>
	auto integrate(UnaryFunction& function, double begin, const double end, const size_t num_steps) {
		const double step = (end - begin) / num_steps;
		auto value = function(begin + 0.5 * step);
		size_t iterNum = 1U;
		do {
			value += function(begin + (iterNum + 0.5) * step);
		} while (++iterNum < num_steps);
		return step * value;
	}

	std::complex<double> complex_complete_elliptic_integral(std::complex<double> k) {
		//if (std::abs(k.imag()) < 1e-12) {
		//	return boost::math::ellint_1(k.real());
		//}

		auto integrand = [k](const double theta) -> std::complex<double> {
			std::complex<double> k_sin = k * sin(theta);
			return (1. / (1. - (k_sin * k_sin)));
		};

		return boost::math::quadrature::gauss_kronrod<double, 61>::integrate(integrand, 0., M_PI_2, 10, 1e-9);
	};

	template <size_t Max_Power>
	std::complex<double> elliptic_power_series(std::complex<double> k) {
		std::complex<double> value{1., 0.};
		k *= k;
		for (size_t n = 1U; n < Max_Power; ++n)
		{
			double coeff = boost::math::double_factorial<double>(2 * n - 1) / boost::math::double_factorial<double>(2 * n);
			value += boost::math::pow<2>(coeff) * std::pow(k, n);
		}
		return M_PI_2 * value;
	}

	void SimpleCubic::computeValues()
	{
		/*
		* Compare
		*/
		constexpr size_t NUM = 100;
		constexpr double DELTA = 1. / NUM;
		std::cout << std::scientific << std::setprecision(15);
		std::cout << "BOOST:\n\n";
		for (size_t i = 0U; i < NUM; ++i)
		{
			std::cout << boost::math::ellint_1(i * DELTA) << ";";
		}
		std::cout << "\n\n\nBoost Integration:\n\n";
		for (size_t i = 0U; i < NUM; ++i)
		{
			std::cout << complex_complete_elliptic_integral(i * DELTA).real() << ";";
		}
		std::cout << "\n\n\nPower series n = 10:\n\n";
		for (size_t i = 0U; i < NUM; ++i)
		{
			std::cout << elliptic_power_series<10>(i * DELTA).real() << ";";
		}
		std::cout << "\n\n\nPower series n = 20:\n\n";
		for (size_t i = 0U; i < NUM; ++i)
		{
			std::cout << elliptic_power_series<20>(i * DELTA).real() << ";";
		}
		std::cout << "\n\n\nPower series n = 100:\n\n";
		for (size_t i = 0U; i < NUM; ++i)
		{
			std::cout << elliptic_power_series<100>(i * DELTA).real() << ";";
		}
		std::cout << std::endl;

		return;
		/*
		*  Algorithm
		*/
		step = 3. / Constants::BASIS_SIZE;
		values.resize(Constants::BASIS_SIZE);
		double gamma = 0.5 * step;
		auto elliptic_integral = [&gamma](const double phi) -> std::complex<double> {
			const std::complex<double> buf = 0.5 * (gamma + cos(phi));
			const std::complex<double> k = sqrt(1. - buf * buf);

			return complex_complete_elliptic_integral(k);
		};

		std::vector<std::complex<double>> returned_values(Constants::BASIS_SIZE);
		for (int g = 0; g < Constants::BASIS_SIZE; ++g)
		{
			gamma = (0.5 + g) * step;
			returned_values[g] = (2. * M_1_PI * M_1_PI * M_1_PI) * integrate(elliptic_integral, -M_PI, M_PI, Constants::BASIS_SIZE);
		}
		computed = true;
	}
}