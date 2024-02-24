#pragma once
#include <string>
#include <vector>
#include "GlobalDefinitions.hpp"
#define TWO_K_DISC (2 * Constants::K_DISCRETIZATION)

namespace Hubbard {
	namespace Constants {
		extern int K_DISCRETIZATION;
		extern int BASIS_SIZE;
		extern int HALF_BASIS;
		extern int QUARTER_BASIS;
		extern int EIGHTH_BASIS;

		extern global_floating_type PI_DIV_DISCRETIZATION;

		const std::vector<std::string> option_list{ "T", "U", "V" };

		inline void setBasis(int setTo) {
			BASIS_SIZE = setTo;
			HALF_BASIS = setTo / 2;
			QUARTER_BASIS = setTo / 4;
			EIGHTH_BASIS = setTo / 8;
		};

		inline void setDiscretization(int setTo)
		{
			K_DISCRETIZATION = setTo;
			PI_DIV_DISCRETIZATION = BASE_PI / setTo;
		};
	}
}