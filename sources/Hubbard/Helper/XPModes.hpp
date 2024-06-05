#pragma once
#include "ModeHelper.hpp"
#include <Utility/Numerics/iEoM/XPResolvent.hpp>

namespace Hubbard::Helper {
	class XPModes : public ModeHelper, protected Utility::Numerics::iEoM::XPResolvent<XPModes, global_floating_type>
	{
		friend struct Utility::Numerics::iEoM::XPResolvent<XPModes, global_floating_type>;
	protected:
		using _parent_algorithm = Utility::Numerics::iEoM::XPResolvent<XPModes, global_floating_type>;

		static constexpr size_t hermitian_size = 7U;
		static constexpr size_t antihermitian_size = 5U;

		const std::array<int, hermitian_size> hermitian_offsets;
		const std::array<int, antihermitian_size> antihermitian_offsets;

		static constexpr std::array<int, 4> cdw_basis_positions{ 2,3,9,10 };

		void fill_M();
		virtual void fillMatrices() override;
		void createStartingStates();
	public:
		XPModes(Utility::InputFileReader& input)
			: ModeHelper(input), _parent_algorithm(this, SQRT_SALT),
			hermitian_offsets{
				0,							Constants::BASIS_SIZE,
				2 * Constants::BASIS_SIZE,	(5 * Constants::BASIS_SIZE) / 2,
				3 * Constants::BASIS_SIZE,	4 * Constants::BASIS_SIZE,
				5 * Constants::BASIS_SIZE
			}, antihermitian_offsets{
				0,									Constants::BASIS_SIZE,
				2 * Constants::BASIS_SIZE,			(5 * Constants::BASIS_SIZE) / 2,
				3 * Constants::BASIS_SIZE
			}
		{ };

		virtual bool matrix_is_negative() override;
		virtual std::vector<ResolventReturnData> computeCollectiveModes() override;
	};
}