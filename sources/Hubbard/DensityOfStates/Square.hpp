#pragma once
#include "BaseDOS.hpp"

namespace Hubbard::DensityOfStates {
	struct Square : public BaseDOS {
		virtual void computeValues() override;
	};
}