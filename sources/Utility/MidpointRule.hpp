#pragma once
#include <utility>

namespace Utility::NumericalSolver::Integration{
    template <class UnaryFunction, class... Args>
	auto midpoint_rule(UnaryFunction& function, double begin, const double end, const size_t num_steps, Args&&... args) {
		const double step = (end - begin) / num_steps;
		auto value = function(begin + 0.5 * step, std::forward<Args>(args)...);
		size_t iterNum = 1U;
		do {
			value += function(begin + (iterNum + 0.5) * step, std::forward<Args>(args)...);
		} while (++iterNum < num_steps);
		return step * value;
	}

	template <class UnaryFunction>
	auto midpoint_rule(UnaryFunction& function, double begin, const double end, const size_t num_steps) {
		const double step = (end - begin) / num_steps;
		auto value = function(begin + 0.5 * step);
		size_t iterNum = 1U;
		do {
			value += function(begin + (iterNum + 0.5) * step);
		} while (++iterNum < num_steps);
		return step * value;
	}
}