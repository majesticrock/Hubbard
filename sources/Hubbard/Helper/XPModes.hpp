#pragma once
#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	class XPModes : public ModeHelper
	{
	protected:
		struct matrix_wrapper {
			Matrix_L eigenvectors;
			Vector_L eigenvalues;

			inline matrix_wrapper() {};

			inline explicit matrix_wrapper(Eigen::Index size)
				: eigenvectors(Matrix_L::Zero(size, size)), eigenvalues(Vector_L::Zero(size))
			{};

			inline Matrix_L reconstruct_matrix() const
			{
				return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.adjoint();
			};

			static matrix_wrapper pivot_and_solve(Matrix_L& toSolve);
			static bool is_non_negative(Matrix_L& toSolve);
		};

		static constexpr size_t hermitian_size = 7U;
		static constexpr size_t antihermitian_size = 5U;

		Matrix_L K_plus, K_minus, L;
		std::array<Vector_L, 2> startingState_SC;
		std::array<Vector_L, 2> startingState_CDW;
		std::array<Vector_L, 2> startingState_AFM;
		std::array<Vector_L, 2> startingState_AFM_transversal;

		const std::array<int, hermitian_size> hermitian_offsets;
		const std::array<int, antihermitian_size> antihermitian_offsets;

		static constexpr std::array<int, 4> cdw_basis_positions{ 2,3,9,10 };

		virtual void fillMatrices() override;
		void createStartingStates();
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