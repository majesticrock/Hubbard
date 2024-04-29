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

		virtual bool matrix_is_negative() override;
		virtual std::vector<ResolventReturnData> computeCollectiveModes() override;
	};
}