#include "BaseDOS.hpp"

namespace Hubbard::DensityOfStates {
	std::vector<double> BaseDOS::values;
	bool BaseDOS::computed = false;
	double BaseDOS::step = 0;
}