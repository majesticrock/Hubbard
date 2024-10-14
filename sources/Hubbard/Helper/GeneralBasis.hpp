#pragma once
#include "ModeHelper.hpp"
#include <Utility/Numerics/iEoM/GeneralResolvent.hpp>

namespace Hubbard::Helper {
	class GeneralBasis : public ModeHelper, protected Utility::Numerics::iEoM::GeneralResolvent<GeneralBasis, global_floating_type>
	{
		friend struct Utility::Numerics::iEoM::GeneralResolvent<GeneralBasis, global_floating_type>;
	protected:
		using _parent_algorithm = Utility::Numerics::iEoM::GeneralResolvent<GeneralBasis, global_floating_type>;

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