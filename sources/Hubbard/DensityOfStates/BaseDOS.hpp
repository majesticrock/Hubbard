#pragma once
#include <vector>
#include <algorithm>
#include <cmath>

namespace Hubbard::DensityOfStates {
	// Computes sqrt(1 - x^2)
	template <class RealType>
	inline RealType sqrt_1_minus_x_squared(RealType x) {
		return std::sqrt(1 - x * x);
	}

	template <bool includesZero, class DataType>
	inline void symmetrizeVector(std::vector<DataType>& v) {
		std::reverse(v.begin(), v.end());
		if constexpr (includesZero) {
			v.insert(v.end(), v.rbegin() + 1, v.rend());
		}
		else {
			v.insert(v.end(), v.rbegin(), v.rend());
		}
	}

	struct BaseDOS {
		// Contains the values for gamma in (-2, 0) as the the positive part is symmetric.
		// Open boundaries, as the dos at 0 contains a singularity for d=2
		static std::vector<double> values;
		static double step;
		static bool computed;

		static void renormalize();
		static void printValues();
		static double getNorm();

		virtual void computeValues() = 0;
		virtual ~BaseDOS() = default;
	};
}