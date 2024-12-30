#pragma once
#include "ModeHelper.hpp"
#include <mrock/utility/Numerics/iEoM/GeneralResolvent.hpp>

namespace Hubbard::Helper {
	class GeneralBasis : public ModeHelper, protected mrock::utility::Numerics::iEoM::GeneralResolvent<GeneralBasis, global_floating_type>
	{
		friend struct mrock::utility::Numerics::iEoM::GeneralResolvent<GeneralBasis, global_floating_type>;
	protected:
		using _parent_algorithm = mrock::utility::Numerics::iEoM::GeneralResolvent<GeneralBasis, global_floating_type>;

		void fill_M();
		virtual void fillMatrices() override;
		void createStartingStates();

		//Debug functions
		void printM(int i, int j) const;
		void printMomentumBlocks() const;
		void printDOSBlocks() const;

		inline Eigen::Index get_index(int basis_term, int momentum_index) const {
			return basis_term * Constants::BASIS_SIZE + momentum_index;
		}
	public:
		GeneralBasis(mrock::utility::InputFileReader& input)
			: ModeHelper(input), _parent_algorithm(this, SQRT_SALT) { };

		virtual bool matrix_is_negative() override;
		virtual std::vector<ResolventReturnData> computeCollectiveModes() override;
	};
}