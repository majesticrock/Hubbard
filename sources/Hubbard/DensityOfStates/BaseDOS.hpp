#pragma once
#include <vector>

namespace Hubbard::DensityOfStates {
	struct BaseDOS {
		// Contains the values for gamma in [-2, 0) as the the positive part is symmetric.
		static std::vector<double> values;
		static double step;
		static bool computed;
		virtual void computeValues() = 0;
	};
}