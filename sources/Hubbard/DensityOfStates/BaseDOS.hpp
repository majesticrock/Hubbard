#pragma once
#include <vector>

namespace Hubbard::DensityOfStates {
	struct BaseDOS {
		// Contains the values for gamma in (-2, 0) as the the positive part is symmetric.
		// Open boundaries, as the dos at 0 contains a singularity for d=2
		static std::vector<double> values;
		static double step;
		static bool computed;
		static void renormalize();
		static void printValues();

		virtual void computeValues() = 0;
		virtual ~BaseDOS() = default;
	};
}