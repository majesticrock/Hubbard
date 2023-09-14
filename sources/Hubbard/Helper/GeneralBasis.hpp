#pragma once
#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	class GeneralBasis : public ModeHelper
	{
	protected:
		MatrixCL M, N;

		virtual void fillMatrices() override;

		//Debug functions
		void printM(int i, int j) const;
		void printMomentumBlocks() const;
		void printDOSBlocks() const;
	public:
		GeneralBasis(Utility::InputFileReader& input) : ModeHelper(input) { };

		virtual std::vector<Resolvent_L> computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) override;
	};
}