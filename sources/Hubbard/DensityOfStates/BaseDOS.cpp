#include "BaseDOS.hpp"
#include <iostream>
#include <numeric>

namespace Hubbard::DensityOfStates {
	std::vector<dos_precision> BaseDOS::values;
	std::vector<abscissa_t> BaseDOS::abscissa;
	bool BaseDOS::computed = false;
	dos_precision BaseDOS::step = 0;

	dos_precision BaseDOS::integrateValues()
	{
		return step * std::reduce(values.begin(), values.end());
	}

	void BaseDOS::printValues()
	{
		for (const auto& val : values) {
			std::cout << val << " ";
		}
		std::cout << std::endl;
	}
}