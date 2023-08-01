#pragma once
#include "BaseDOS.hpp"

namespace Hubbard::DensityOfStates {
	struct SimpleCubic : public BaseDOS {
		static constexpr double LOWER_BORDER = -3;
		virtual void computeValues() override;

		template <class ResultType>
		class DOSIntegrator {
		private:
			ResultType result;
			ResultType buffer;

			template <bool byValue, class UnaryFunction>
			const ResultType& _internal_integrate(const UnaryFunction& F) {
				if constexpr (byValue) {
					result = values[0] * F(LOWER_BORDER);
				}
				else {
					F(LOWER_BORDER, buffer);
					result = values[0] * buffer;
				}
				for (size_t i = 1U; i < values.size(); ++i)
				{
					if constexpr (byValue) {
						result += values[i] * F(step * i + LOWER_BORDER);
					}
					else {
						F(step * i + LOWER_BORDER, buffer);
						result += values[i] * buffer;
					}
				}
				result *= step;
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
	};
}