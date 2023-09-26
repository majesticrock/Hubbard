#pragma once
#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	class XPModes : public ModeHelper
	{
	protected:
		Matrix_L K_plus, K_minus, L;
		virtual void fillMatrices() override;
	public:
		XPModes(Utility::InputFileReader& input) : ModeHelper(input) { };

		virtual std::vector<ResolventReturnData> computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) override;
	};
}