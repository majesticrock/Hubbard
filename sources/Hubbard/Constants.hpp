#pragma once
#include <string>
#include <vector>
#include "GlobalDefinitions.hpp"

namespace Hubbard {
	namespace Constants {
		extern int K_DISCRETIZATION;
		extern int BASIS_SIZE;
		extern int HALF_BASIS;
		extern global_floating_type PI_DIV_DISCRETIZATION;

		const std::vector<std::string> option_list{ "T", "U", "V" };

		inline void setDiscretization(int setTo) 
		{
			K_DISCRETIZATION = setTo;
			BASIS_SIZE = 4 * setTo * setTo;
			HALF_BASIS = BASIS_SIZE / 2;
			PI_DIV_DISCRETIZATION = BASE_PI / setTo;
		};
	}
}