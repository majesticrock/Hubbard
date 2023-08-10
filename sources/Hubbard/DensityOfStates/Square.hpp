#pragma once
#include <boost/multiprecision/float128.hpp>
#include "BaseDOS.hpp"
#include <cmath>

namespace Hubbard::DensityOfStates {
	typedef boost::multiprecision::float128 abscissa_t;
	struct Square : public BaseDOS {

		static std::vector<double> regular_values;
		static std::vector<double> singular_values_linear;
		static std::vector<double> singular_values_quadratic;
		static std::vector<double> singular_weights;

		static std::vector<abscissa_t> abscissa;
		static std::vector<abscissa_t> upper_border_to_abscissa;
		static std::vector<double> weights;

		static double LOWER_BORDER;
		static double b_minus_a_halved;

		static void tanh_sinh();

		template <class ResultType>
		static inline ResultType slope_m(const ResultType& previous_function_value, const ResultType& current_function_value) {
			return (current_function_value - previous_function_value) / step;
		};
		static inline double neighbouringDifference_linear(size_t i) {
			return (singular_values_linear[i] - singular_values_linear[i - 1]);
		};
		static inline double neighbouringDifference_quadratic(size_t i) {
			return (singular_values_quadratic[i] - singular_values_quadratic[i - 1]);
		};
		virtual void computeValues() override;
		
		template <class ResultType>
		class DOSIntegrator {
		private:
			ResultType result;
			ResultType buffer;

			template <bool byValue, class UnaryFunction>
			const ResultType& _internal_integrate(const UnaryFunction& F) {
				if constexpr (byValue) {
					result = values[0] * weights[0] * F(static_cast<double>(abscissa.front()));
				}
				else {
					F(static_cast<double>(abscissa.front()), buffer);
					result = values[0] * weights[0] * buffer;
				}
				for (size_t i = 1U; i < values.size(); ++i)
				{
					if constexpr (byValue) {
						result += values[i] * weights[i] * F(static_cast<double>(abscissa[i]));
					}
					else {
						F(static_cast<double>(abscissa[i]), buffer);
						result += values[i] * weights[i] * buffer;
					}
				}
				result *= 2;
				return result;
			};
		public:
			// This function passes the result of F by reference,
			// i.e. F(gamma, result) and expects F to fill result accordingly.
			template <class UnaryFunction>
			inline const ResultType& integrate_by_reference(const UnaryFunction& F) {
				return _internal_integrate<false, UnaryFunction>(F);
			};

			// This function assumes that F returns its result by value.
			template <class UnaryFunction>
			inline const ResultType& integrate_by_value(const UnaryFunction& F) {
				return _internal_integrate<true, UnaryFunction>(F);
			};

			DOSIntegrator() = default;
			// If it is neccessary for the ResultType to be initaliazed,
			// e.g. give a vector a certain size. ResultType needs to have a copy constructor though
			DOSIntegrator(const ResultType& initialize_result)
				: result(initialize_result), buffer(initialize_result)
			{};
		};
		/*class DOSIntegrator {
		private:
			ResultType previous_function_value; // F(a_n)
			ResultType current_function_value; // F(a_n+1)
			ResultType result;
			ResultType singular_value;

			double gamma{};

			template <bool byValue, class UnaryFunction>
			const ResultType& _internal_integrate(const UnaryFunction& F) {
				gamma = LOWER_BORDER + step;
				if constexpr (byValue) {
					previous_function_value = F(LOWER_BORDER); // F(a_n)
					current_function_value = F(gamma); // F(a_n+1)
				}
				else {
					F(LOWER_BORDER, previous_function_value); // F(a_n)
					F(gamma, current_function_value); // F(a_n+1)
				}
				result = regular_values[0] * previous_function_value + regular_values[1] * current_function_value;

				//singular_value = neighbouringDifference_linear(1) * (previous_function_value - LOWER_BORDER * slope_m(previous_function_value, current_function_value))
				//	+ neighbouringDifference_quadratic(1) * slope_m(previous_function_value, current_function_value);
				singular_value = singular_weights[0] * previous_function_value + singular_weights[1] * current_function_value;

				for (size_t i = 2U; i < regular_values.size(); ++i)
				{
					previous_function_value = current_function_value;
					gamma += step;
					if constexpr (byValue) {
						current_function_value = F(gamma); // F(a_n+1)
					}
					else {
						F(gamma, current_function_value); // F(a_n+1)
					}

					result += regular_values[i] * current_function_value; // Integral R(alpha) * F(gamma)

					singular_value += singular_weights[i] * current_function_value;
					//singular_value += neighbouringDifference_linear(i)
					//	* (previous_function_value - (gamma - step) * slope_m(previous_function_value, current_function_value));
					//singular_value += neighbouringDifference_quadratic(i) * slope_m(previous_function_value, current_function_value);
				}

				result *= step;
				result -= 0.5 * singular_value;
				return result;
			}
		public:
			// This function passes the result of F by reference,
			// i.e. F(gamma, result) and expects F to fill result accordingly.
			template <class UnaryFunction>
			inline const ResultType& integrate_by_reference(const UnaryFunction& F) {
				return _internal_integrate<false, UnaryFunction>(F);
			};

			// This function assumes that F returns its result by value.
			template <class UnaryFunction>
			inline const ResultType& integrate_by_value(const UnaryFunction& F) {
				return _internal_integrate<true, UnaryFunction>(F);
			};

			DOSIntegrator() = default;
			// If it is neccessary for the ResultType to be initaliazed,
			// e.g. give a vector a certain size. ResultType needs to have a copy constructor though
			DOSIntegrator(const ResultType& initialize_result) 
				: previous_function_value(initialize_result), current_function_value(initialize_result), 
				result(initialize_result), singular_value(initialize_result)
			{};
		};*/
	};
}