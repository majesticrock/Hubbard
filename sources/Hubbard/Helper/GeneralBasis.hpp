#pragma once
#include "ModeHelper.hpp"
#include <Utility/Numerics/iEoM/GeneralResolvent.hpp>

namespace Hubbard::Helper {
	class GeneralBasis : public ModeHelper, protected Utility::Numerics::iEoM::GeneralResolvent<GeneralBasis, complex_prec>
	{
		friend struct Utility::Numerics::iEoM::GeneralResolvent<GeneralBasis, complex_prec>;
	protected:
		using _parent_algorithm = Utility::Numerics::iEoM::GeneralResolvent<GeneralBasis, complex_prec>;

		static constexpr int NUMBER_OF_GREENSFUNCTIONS = 4;

		void fill_M();
		virtual void fillMatrices() override;
		void createStartingStates();

		//Debug functions
		void printM(int i, int j) const;
		void printMomentumBlocks() const;
		void printDOSBlocks() const;
	public:
		GeneralBasis(Utility::InputFileReader& input)
			: ModeHelper(input), _parent_algorithm(this, SQRT_SALT) { };

		virtual bool matrix_is_negative() override;
		virtual std::vector<ResolventReturnData> computeCollectiveModes() override;
	};
}