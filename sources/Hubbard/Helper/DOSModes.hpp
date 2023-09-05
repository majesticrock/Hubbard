#pragma once
#include "GeneralBasis.hpp"

namespace Hubbard::Helper {
	class DOSModes : public GeneralBasis
	{
	private:
		virtual void fillBlock(int i, int j) override;

	public:
		DOSModes(Utility::InputFileReader& input) : GeneralBasis(input) {};
	};
}