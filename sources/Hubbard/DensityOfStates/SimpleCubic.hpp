#pragma once
#include "BaseDOS.hpp"

namespace Hubbard::DensityOfStates {
	class SimpleCubic : public BaseDOS {
	public:
		virtual void computeValues() override;
	};
}