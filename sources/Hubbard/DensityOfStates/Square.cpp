#include "Square.hpp"
#include "../Constants.hpp"
#include <numeric>
#include <boost/math/constants/constants.hpp>

using boost::math::constants::pi;
using boost::math::constants::half_pi;
using boost::multiprecision::exp;
using boost::multiprecision::sinh;
using boost::multiprecision::cosh;
using boost::multiprecision::tanh;
using boost::multiprecision::pow;

namespace Hubbard::DensityOfStates {
	std::vector<abscissa_t> Square::abscissa;
	std::vector<abscissa_t> Square::upper_border_to_abscissa;
	std::vector<double> Square::weights;

	std::vector<double> Square::regular_values;
	std::vector<double> Square::singular_values_linear;
	std::vector<double> Square::singular_values_quadratic;

	std::vector<double> Square::singular_weights;
	double Square::LOWER_BORDER;
	double Square::b_minus_a_halved;

	template <class RealType>
	// Returns x*ln((x/2)^2)
	inline RealType x_log_x_2_squared(RealType x) {
		if (std::abs(x) < CUT_OFF) return 0;
		return x * (std::log(x * x) - LOG_4);
	}

	template <class T>
	void expand_vector(std::vector<T>& vec) {
		if (vec.size() == 0) {
			vec.push_back(T{});
			return;
		}
		int n = vec.size();
		vec.reserve(2 * n - 1);
		for (int i = n - 1; i > 0; --i)
		{
			vec.insert(vec.begin() + i, T{});
		}
	}

	template <class RealType>
	struct tanh_sinh_helper {
	private:
		RealType _sinh_x{};
		// Integrating from 0 to 2
		const double _lower_border{ 0 };
		const double _upper_border{ 2 };
		const double _center{ (_upper_border + _lower_border) / 2 }; // (a+b)/2
		const double _half_distance{ (_upper_border - _lower_border) / 2 }; // (b-a)/2
	public:
		inline void set_sinh_x(RealType x) {
			_sinh_x = sinh(x);
		};
		inline RealType compute_abscissa() const {
			return _upper_border - compute_upper_border_to_abscissa();// _center + _half_distance * std::tanh(LONG_PI_2 * _sinh_x);
		};
		// Computes b - x
		inline RealType compute_upper_border_to_abscissa() const {
			return _half_distance * 2. / 
				(1 + exp(pi<abscissa_t>() * _sinh_x));
		};
		// Computes a - x
		inline RealType compute_lower_border_to_abscissa() const {
			return -2. * _half_distance * exp(pi<abscissa_t>() * _sinh_x) / (1 + exp(pi<abscissa_t>() * _sinh_x));
		};
		inline double compute_weight(int k) const {
			return (Square::step * LONG_PI_2 * std::cosh(k * Square::step)) 
				/ static_cast<long double>(pow(cosh(half_pi<abscissa_t>() * _sinh_x), 2));
		};

		inline decltype(_half_distance) half_distance() const {
			return _half_distance;
		};

		tanh_sinh_helper() = default;
		tanh_sinh_helper(RealType a, RealType b)
			: _lower_border(a), _upper_border(b), _center((b + a) / 2), _half_distance((b - a) / 2) {};
	};

	void Square::tanh_sinh() {
		auto compute_DOS = [](abscissa_t gamma, abscissa_t one_minus_gamma) -> double {
			gamma *= 0.5;
			one_minus_gamma *= 0.5;
			return (LONG_1_PI * LONG_1_PI)
				* boost::math::ellint_1((gamma < 0.25 ? sqrt(1 - gamma * gamma) : sqrt(one_minus_gamma * (1 + gamma)))).convert_to<double>();
		};

		tanh_sinh_helper<abscissa_t> tsh{ 0, 2 };
		double old_integral{};

		auto fill_vectors = [&](int k, bool sign) {
			do {
				tsh.set_sinh_x((sign ? k : -k) * step);
				upper_border_to_abscissa.push_back(tsh.compute_upper_border_to_abscissa());
				abscissa.push_back(tsh.compute_abscissa());
				weights.push_back(tsh.compute_weight((sign ? k : -k)));
				values.push_back(compute_DOS(abscissa.back(), upper_border_to_abscissa.back()));
				++k;

				old_integral += values.back() * weights.back();
			} while (std::abs(values.back() * weights.back()) > 1e-11 || std::abs(weights.back()) > 1e-8);
			return k - 1;
		};

		int min_k{ -fill_vectors(0, false) };

		std::reverse(values.begin(), values.end());
		std::reverse(weights.begin(), weights.end());
		std::reverse(abscissa.begin(), abscissa.end());
		std::reverse(upper_border_to_abscissa.begin(), upper_border_to_abscissa.end());

		fill_vectors(1, true);
		old_integral *= tsh.half_distance();

		double new_integral{};
		double error{ 100.0 };
		size_t level{};
		while (error > 1e-12) {
			++level;
			step /= 2;
			expand_vector(values);
			expand_vector(abscissa);
			expand_vector(upper_border_to_abscissa);
			expand_vector(weights);
			min_k *= 2;

			new_integral = 0;
			for (int k = 0; k < values.size(); ++k)
			{
				if (k % 2 != 0) {
					tsh.set_sinh_x((k + min_k) * step);
					abscissa[k] = tsh.compute_abscissa();
					upper_border_to_abscissa[k] = tsh.compute_upper_border_to_abscissa();
					weights[k] = tsh.compute_weight(k + min_k);
					values[k] = compute_DOS(abscissa[k], upper_border_to_abscissa[k]);
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

		std::cout << "Exit after " << level << " levels with error = " << std::abs(0.5 - new_integral) << std::endl;
		std::cout << "Total amount of values = " << values.size() << std::endl;
	}

	void Square::computeValues()
	{
		step = std::ldexp(1, -1);
		tanh_sinh();

		//step = 2. / Constants::BASIS_SIZE;
		//LOWER_BORDER = -2 + 0.5 * step;
		//
		//values.reserve(2 * Constants::BASIS_SIZE);
		//values.resize(Constants::BASIS_SIZE);
		//
		//for (int g = 0; g < Constants::BASIS_SIZE; ++g)
		//{
		//	const long double gamma = (0.5 + g) * step;
		//	values[g] = (LONG_1_PI * LONG_1_PI)
		//		* boost::math::ellint_1(sqrt_1_minus_x_squared(0.5 * gamma));
		//}
		//symmetrizeVector<false>(values);

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

		auto weight = [](double gamma) {
			if (gamma > -1e-12 && gamma < 1e-12) return 0.0;
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