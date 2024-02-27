#pragma once
#include <Eigen/Dense>

namespace Utility {
	// Pivots a matrix so that all offdiagonal 0 blocks are contiguous
	// The permutation matrix is returned
	// epsilon is used to determine if a matrix element is 0 (especially important for floating point operations)
	template<class EigenMatrixType>
	auto pivot_to_block_structure(const EigenMatrixType& matrix, const typename EigenMatrixType::Scalar epsilon = 1e-12) {
		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P(matrix.rows());
		P.setIdentity();

		auto permuted_matrix_element = [&matrix, &P](Eigen::Index i, Eigen::Index j) {
			return matrix.coeffRef(P.indices()(i), P.indices()(j));
			};

		for (int i = 0; i < matrix.rows(); ++i) {
			int offset = 1;
			for (int j = i + 1; j < matrix.cols(); ++j) {
				if (permuted_matrix_element(i, j) > epsilon) {
					P.applyTranspositionOnTheRight(j, i + offset);
					++offset;
				}
			}
		}
		return P;
	};

	struct HermitianBlock{
		Eigen::Index position{};
		Eigen::Index size{};
	};

	inline std::ostream& operator<<(std::ostream& os, const HermitianBlock& block){
		os << block.position << "\t" << block.size;
		return os;
	};

	template<class EigenMatrixType>
	auto identify_hermitian_blocks(const EigenMatrixType& matrix, const typename EigenMatrixType::Scalar epsilon = 1e-12) {
		Eigen::Index block_index{ 1 };
		Eigen::Index block_size{};
		std::vector<HermitianBlock> block_indices;
		for (Eigen::Index i = 0; i < matrix.rows(); ++i)
		{
			for (Eigen::Index j = matrix.cols() - 1; j > block_index; --j)
			{
				if (abs(matrix(i, j)) > 1e-12) {
					block_index = j;
					break;
				}
			}
			if (block_index == i) {
				if (block_indices.empty()) {
					block_indices.push_back({ Eigen::Index{}, block_index + 1 });
				}
				else {
					block_size = block_index - (block_indices.back().size + block_indices.back().position) + 1;
					block_indices.push_back({ block_index - block_size + 1, block_size });
				}
			}
		}
		return block_indices;
	};
}