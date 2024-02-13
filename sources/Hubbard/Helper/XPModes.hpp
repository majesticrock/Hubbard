#pragma once
#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	class XPModes : public ModeHelper
	{
	protected:
		Matrix_L K_plus, K_minus, L;
		const std::array<int, 6> hermitian_offsets;
		const std::array<int, 4> antihermitian_offsets;

		static constexpr std::array<int, 4> cdw_basis_positions{ 2,3,8,9 };

		virtual void fillMatrices() override;
	public:
		XPModes(Utility::InputFileReader& input) : ModeHelper(input),
			hermitian_offsets{
				0,							Constants::BASIS_SIZE,
				2 * Constants::BASIS_SIZE,	(5 * Constants::BASIS_SIZE) / 2,
				3 * Constants::BASIS_SIZE,	4 * Constants::BASIS_SIZE
		}, antihermitian_offsets{
			0,									Constants::BASIS_SIZE,
			2 * Constants::BASIS_SIZE,			(5 * Constants::BASIS_SIZE) / 2
		}
		{ };

		virtual bool matrix_is_negative() override;
		virtual std::vector<ResolventReturnData> computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) override;
	};
}