#pragma once
#include <utility>

namespace Utility::NumericalSolver::Integration {
	template <class UnaryFunction, class... Args>
	auto midpoint_rule(UnaryFunction& function, double begin, const double end, unsigned long num_steps, Args&&... args) {
		const double step = (end - begin) / num_steps;
		begin += 0.5 * step;
		auto value = function(begin, std::forward<Args>(args)...);
		while (--num_steps > 0)
		{
			begin += step;
			value += function(begin, std::forward<Args>(args)...);
		}
		return step * value;
	}

	template <class UnaryFunction>
	auto midpoint_rule(UnaryFunction& function, double begin, const double end, unsigned long num_steps) {
		const double step = (end - begin) / num_steps;
		begin += 0.5 * step;
		auto value = function(begin);
		while (--num_steps > 0) {
			begin += step;
			value += function(begin);
		}
		return step * value;
	}
}