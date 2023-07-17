#include "BaseDOS.hpp"
#include <iostream>

namespace Hubbard::DensityOfStates {
	std::vector<double> BaseDOS::values;
	bool BaseDOS::computed = false;
	double BaseDOS::step = 0;
	void BaseDOS::printValues()
	{
		for (const auto& val : values) {
			std::cout << val << " ";
		}
		std::cout << std::endl;
	}
}