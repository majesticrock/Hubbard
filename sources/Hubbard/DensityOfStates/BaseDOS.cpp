#include "BaseDOS.hpp"
#include <iostream>
#include <numeric>

namespace Hubbard::DensityOfStates {
	std::vector<dos_precision> BaseDOS::values;
	std::vector<abscissa_t> BaseDOS::abscissa;
	std::vector<dos_precision> BaseDOS::weights;
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
	void BaseDOS::clearAll()
	{
		values.clear();
		abscissa.clear();
		weights.clear();
		computed = false;
		step = 0;
	}
	void BaseDOS::printValuesAndAbscissa()
	{
		for (size_t i = 0U; i < values.size(); ++i)
		{
			std::cout << abscissa[i] << "; " << values[i] << std::endl;
		}
		std::cout << std::endl;
	}
}