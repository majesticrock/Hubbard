#pragma once
#include "BaseDOS.hpp"

namespace Hubbard::DensityOfStates {
	struct SimpleCubic : public BaseDOS {
	private:
		static constexpr double LOWER_BORDER = -3;
	public:
		template <class UnaryFunction, class ResultType>
		static inline void integrate_with_dos(const UnaryFunction& F, ResultType& result) {
			for (size_t i = 0U; i < values.size(); ++i)
			{
				F(result, step * i + LOWER_BORDER, values[i]);
			}
			result *= step;
		};

		template <class UnaryFunction>
		static inline auto integrate_with_dos(const UnaryFunction& F) {
			decltype(F(double{})) result{values[0] * F(LOWER_BORDER)};

			for (size_t i = 1U; i < values.size(); ++i)
			{
				result += values[i] * F(step * i + LOWER_BORDER);
			}
			result *= step;
			return result;
		};

		virtual void computeValues() override;
	};
}