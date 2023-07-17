#pragma once
#include "BaseDOS.hpp"

namespace Hubbard::DensityOfStates {
	class SimpleCubic : public BaseDOS {
		virtual void computeValues() override;
	};
}