#pragma once
#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	class XPModes : public ModeHelper
	{
	protected:
		static constexpr size_t hermitian_size = 7U;
		static constexpr size_t antihermitian_size = 5U;

		Matrix_L K_plus, K_minus, L;
		const std::array<int, hermitian_size> hermitian_offsets;
		const std::array<int, antihermitian_size> antihermitian_offsets;

		static constexpr std::array<int, 4> cdw_basis_positions{ 2,3,9,10 };

		virtual void fillMatrices() override;
	public:
		XPModes(Utility::InputFileReader& input) : ModeHelper(input),
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
		virtual std::vector<ResolventReturnData> computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) override;
	};
}