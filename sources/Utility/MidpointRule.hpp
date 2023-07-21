#pragma once
#include <utility>

namespace Utility::NumericalSolver::Integration {
	template <class UnaryFunction, class RealType, class... Args>
	auto midpoint_rule(UnaryFunction& function, RealType begin, const RealType end, unsigned long num_steps, Args&&... args) {
		const RealType step = (end - begin) / num_steps;
		begin += 0.5 * step;
		auto value = function(begin, std::forward<Args>(args)...);
		while (--num_steps > 0)
		{
			begin += step;
			value += function(begin, std::forward<Args>(args)...);
		}
		return step * value;
	}

	template <class UnaryFunction, class RealType>
	auto midpoint_rule(UnaryFunction& function, RealType begin, const RealType end, unsigned long num_steps) {
		const RealType step = (end - begin) / num_steps;
		begin += 0.5 * step;
		auto value = function(begin);
		while (--num_steps > 0) {
			begin += step;
			value += function(begin);
		}
		return step * value;
	}
}