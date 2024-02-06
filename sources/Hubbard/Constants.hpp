#pragma once
#include <string>
#include <vector>
#include "GlobalDefinitions.hpp"

namespace Hubbard {
	namespace Constants {
		extern int K_DISCRETIZATION;
		extern int BASIS_SIZE;
		extern int HALF_BASIS;
		extern size_t SPINOR_SIZE;

		extern global_floating_type PI_DIV_DISCRETIZATION;

		const std::vector<std::string> option_list{ "T", "U", "V" };

		inline void setBasis(int setTo) {
			BASIS_SIZE = setTo;
			HALF_BASIS = setTo / 2;
		};

		inline void setDiscretization(int setTo)
		{
			K_DISCRETIZATION = setTo;
			PI_DIV_DISCRETIZATION = BASE_PI / setTo;
		};
	}
}