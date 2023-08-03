#include "Square.hpp"
#include "../Constants.hpp"
#include <numeric>

namespace Hubbard::DensityOfStates {
	std::vector<double> Square::abscissa;
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

	void expand_vector(std::vector<double>& vec){
		if(vec.size() == 0){
			vec.push_back(0);
			return;
		}
		int n = vec.size();
		vec.reserve(2 * n - 1);
		for (int i = n - 1; i > 0; --i)
		{
			vec.insert(vec.begin() + i, 0);
		}
	}

	void Square::tanh_sinh(){
		auto compute_DOS = [](long double gamma){
			//std::cout << std::scientific << std::setprecision(12) << "         t = " << sqrt_1_minus_x_squared<long double>(0.5 * gamma) << std::endl;
			return (LONG_1_PI * LONG_1_PI) * boost::math::ellint_1(sqrt_1_minus_x_squared(0.5L * gamma));
		};
		auto compute_abscissa = [](int k) {
			return 0.5 + 0.5 * std::tanh(0.5 * LONG_PI * std::sinh(k * step));
		};
		auto compute_weight = [](int k){
			return 0.5 * 0.5 * step * LONG_PI * std::cosh(k * step) / std::pow(std::cosh(0.5 * LONG_PI * std::sinh(k * step)), 2);
		};

		double old_integral = 0;
		int min_k = 0;
		do {
			abscissa.push_back(compute_abscissa(min_k));
			weights.push_back(compute_weight(min_k));
			values.push_back(compute_DOS(abscissa.back()));
			--min_k;

			old_integral += values.back() * weights.back();
		} while (std::abs(values.back() * weights.back()) > 1e-9);
		++min_k;
		std::reverse(values.begin(), values.end());
		int max_k = 1;
		do {
			abscissa.push_back(compute_abscissa(max_k));
			weights.push_back(compute_weight(max_k));

			values.push_back(compute_DOS(abscissa.back()));
			++max_k;

			old_integral += values.back() * weights.back();
		} while (std::abs(values.back() * weights.back()) > 1e-9);

		double new_integral;
		double error = 100;
		size_t level = 0;
		while(error > 1e-10){
			++level;
			step /= 2;
			expand_vector(values);
			expand_vector(abscissa);
			expand_vector(weights);
			min_k *= 2;

			new_integral = 0;
			for (int k = 0; k < values.size(); ++k)
			{
				if(k % 2 != 0){
					abscissa[k] = compute_abscissa(k + min_k);
					weights[k] = compute_weight(k + min_k);
					values[k] = compute_DOS(abscissa[k]);
				}else{
					weights[k] *= 0.5;
				}
				new_integral += values[k] * weights[k];
			}
			
			auto w_it = weights.begin();
			auto a_it = abscissa.begin();
			for(auto v_it = values.begin(); v_it != values.end();){
				if( std::abs((*v_it) * (*w_it)) < 1e-9){
					v_it = values.erase(v_it);
					a_it = abscissa.erase(a_it);
					w_it = weights.erase(w_it);
				}else{
					++v_it;
					++a_it;
					++w_it;
				}
			}

			error = std::abs((new_integral - old_integral));
			std::cout << "New value = " << new_integral << "         Error = " << error << std::endl;
			std::cout << "Total amount of values = " << values.size() << std::endl;
			old_integral = new_integral;
		}

		std::cout << "Exit after " << level << " levels with error = " << error << std::endl;
		std::cout << "Total amount of values = " << values.size() << std::endl;
	}

	void Square::computeValues()
	{
		step = 0.01;
		tanh_sinh();

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