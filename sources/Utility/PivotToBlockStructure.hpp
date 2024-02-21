#pragma once
#include <Eigen/Dense>

namespace Utility{

    // Pivots a matrix so that all offdiagonal 0 blocks are contiguous
    // The matrix is changed in place and the permutation matrix is returned
    template<class EigenMatrixType>
    auto pivot_to_block_structure(EigenMatrixType& matrix, double epsilon=1e-12){
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permMatrix(matrix.rows());
        permMatrix.setIdentity();

        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = i + 1; j < matrix.cols(); ++j) {
                if (abs(matrix(i, j)) > epsilon) {
                    permMatrix.applyTranspositionOnTheRight(j, i + 1);
                    matrix.col(j).swap(matrix.col(i + 1));
                    matrix.row(j).swap(matrix.row(i + 1));
                }
            }
        }
        return permMatrix;
    };

    template<class EigenMatrixType>
    auto identify_hermitian_blocks(const EigenMatrixType& matrix, double epsilon=1e-12){
        Eigen::Index block_index{ 1 };
        Eigen::Index block_size{};
        std::vector<std::pair<Eigen::Index, Eigen::Index>> block_indices;
        for (Eigen::Index i = 0; i < matrix.rows(); ++i)
        {
            for (Eigen::Index j = matrix.cols() - 1; j > block_index; --j)
            {
                if(abs(matrix(i, j)) > 1e-12){
                    block_index = j;
                    break;
                }
            }
            if (block_index == i) {
                if(block_indices.empty()){
                    block_indices.push_back(std::make_pair(Eigen::Index{}, block_index + 1));
                } else {
                    block_size = block_index - (block_indices.back().second + block_indices.back().first) + 1;
                    block_indices.push_back(std::make_pair(block_index - block_size + 1, block_size));
                }
            }
        }
        return block_indices;
    };
}