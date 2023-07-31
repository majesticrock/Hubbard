#pragma once
#include "BaseDOS.hpp"
#include <cmath>

namespace Hubbard::DensityOfStates {
	struct Square : public BaseDOS {
	private:
		static double LOWER_BORDER;
		static std::vector<double> regular_values;
		static std::vector<double> singular_values_linear;
		static std::vector<double> singular_values_quadratic;

		template <class ResultType>
		static inline ResultType slope_m(const ResultType& previous_function_value, const ResultType& current_function_value) {
			return (current_function_value - previous_function_value) / step;
		}
	public:
		template <class UnaryFunction, class ResultType>
		static void integrate_with_dos(const UnaryFunction& F, ResultType& result) {
			for (size_t i = 0U; i < values.size(); ++i)
			{
				F(result, step * i + LOWER_BORDER, values[i]);
			}
			result *= step;
		};

		template <class UnaryFunction>
		static inline auto integrate_with_dos(const UnaryFunction& F) {
			using ReturnType = decltype(F(double{}));
			double gamma{ LOWER_BORDER };

			ReturnType previous_function_value{ F(LOWER_BORDER) }; // F(a_n)
			ReturnType current_function_value{ F(gamma) }; // F(a_n+1)
			ReturnType result{ regular_values[0] * previous_function_value
				+ regular_values[1] * current_function_value };

			ReturnType singular_value{
				(singular_values_linear[1] - singular_values_linear[0])
					* (previous_function_value - LOWER_BORDER * slope_m(previous_function_value, current_function_value))
				+ (singular_values_quadratic[1] - singular_values_quadratic[0])
					* slope_m(previous_function_value, current_function_value) };

			for (size_t i = 2U; i < regular_values.size(); ++i)
			{
				previous_function_value = current_function_value;
				gamma += step;
				current_function_value = F(gamma); // F(a_n+1)

				result += regular_values[i] * current_function_value; // Integral R(alpha) * F(gamma)

				singular_value += (singular_values_linear[i] - singular_values_linear[i - 1])
					* (previous_function_value - (gamma - step) * slope_m(previous_function_value, current_function_value));

				singular_value += (singular_values_quadratic[i] - singular_values_quadratic[i - 1])
					* slope_m(previous_function_value, current_function_value);
			}

			result *= step;
			result -= 0.5 * singular_value;
			return result;
		};

		virtual void computeValues() override;
	};
}