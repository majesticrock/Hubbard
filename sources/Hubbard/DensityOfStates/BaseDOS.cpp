#include "BaseDOS.hpp"
#include <iostream>
#include <numeric>

namespace Hubbard::DensityOfStates {
	std::vector<double> BaseDOS::values;
	bool BaseDOS::computed = false;
	double BaseDOS::step = 0;
	void BaseDOS::renormalize()
	{
		// Factor of 2 because we only computed half the DOS due to symmetry
		double norm = 2 * step * std::accumulate(values.begin(), values.end(), double{});
		std::cout << "DOS-Norm: " << norm << std::endl;
		for (auto& val : values) {
			val /= norm;
		}
	}
	void BaseDOS::printValues()
	{
		for (const auto& val : values) {
			std::cout << val << " ";
		}
		std::cout << std::endl;
	}
}