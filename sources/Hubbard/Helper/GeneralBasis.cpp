#include "GeneralBasis.hpp"
#include <chrono>
#include "../DensityOfStates/Square.hpp"

using DOS = Hubbard::DensityOfStates::Square;

namespace Hubbard::Helper {
	void GeneralBasis::fillMatrices()
	{
		M = MatrixCL::Zero(TOTAL_BASIS, TOTAL_BASIS);
		N = MatrixCL::Zero(TOTAL_BASIS, TOTAL_BASIS);

		//#pragma omp parallel for
		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				fillBlock(i, j);
			}
		}
	}

	std::vector<Resolvent_L> GeneralBasis::computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) {
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		std::chrono::time_point end = std::chrono::steady_clock::now();

		std::cout << std::resetiosflags(std::cout.flags());
		fillMatrices();
		std::cout << std::endl;
		//int off = Constants::BASIS_SIZE ;
		//for (size_t k = 0U; k < off; ++k)
		//{
		//	for (size_t i = 0U; i < 2; ++i)
		//	{
		//		for (size_t j = 0U; j < 2; ++j)
		//		{
		//			std::cout << M(j * off + k, i * off + k).real() << "\t";
		//		}
		//		std::cout << std::endl;
		//	}
		//	std::cout << std::endl;
		//}

		if ((M - M.adjoint()).norm() > 1e-12) {
			throw std::runtime_error("M is not Hermitian!");
		}
		if ((N - N.adjoint()).norm() > 1e-12) {
			throw std::runtime_error("N is not Hermitian!");
		}

		end = std::chrono::steady_clock::now();
		std::cout << "Time for filling of M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		Eigen::SelfAdjointEigenSolver<MatrixCL> M_solver(M);
		Vector_L& evs_M = const_cast<Vector_L&>(M_solver.eigenvalues());
		applyMatrixOperation<OPERATION_NONE>(evs_M);

		auto bufferMatrix = N * M_solver.eigenvectors();
		MatrixCL n_hacek = bufferMatrix
			* evs_M.unaryExpr([](global_floating_type x) { return abs(x < SALT) ? 0 : 1. / x; }).asDiagonal()
			* bufferMatrix.adjoint();

		Eigen::SelfAdjointEigenSolver<MatrixCL> norm_solver(n_hacek);
		Vector_L& evs_norm = const_cast<Vector_L&>(norm_solver.eigenvalues());
		applyMatrixOperation<OPERATION_SQRT>(evs_norm);

		n_hacek = norm_solver.eigenvectors()
			* evs_norm.unaryExpr([](global_floating_type x) { return abs(x < SALT) ? 0 : 1. / x; }).asDiagonal()
			* norm_solver.eigenvectors().adjoint();
		// Starting here M is the adjusted solver matrix (s s hackem)
		M = n_hacek * M_solver.eigenvectors() * evs_M.asDiagonal() * M_solver.eigenvectors().adjoint() * n_hacek;
		// Starting here N is the extra matrix that defines |a> (n s hackem N)
		N.applyOnTheLeft(n_hacek);
		// Starting here h_hacek is its own inverse (defining |b>)
		n_hacek = norm_solver.eigenvectors() * evs_norm.asDiagonal() * norm_solver.eigenvectors().adjoint();

		end = std::chrono::steady_clock::now();
		std::cout << "Time for adjusting of the matrices: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		constexpr int NUMBER_OF_GREENSFUNCTIONS = 4;
		/*
		* 0 - SC Higgs
		* 1 - SC Phase
		* 2 - CDW
		* 3 - AFM
		*/
		std::vector<VectorCL> psis(NUMBER_OF_GREENSFUNCTIONS, VectorCL::Zero(TOTAL_BASIS));
		for (int i = 0; i < Constants::BASIS_SIZE; i++)
		{
			psis[0](i) = 1;
			psis[0](Constants::BASIS_SIZE + i) = 1;

			psis[1](i) = 1;
			psis[1](Constants::BASIS_SIZE + i) = -1;

			if (number_of_basis_terms >= 6) {
				psis[2](4 * Constants::BASIS_SIZE + i) = 1;
				psis[2](5 * Constants::BASIS_SIZE + i) = 1;

				psis[3](4 * Constants::BASIS_SIZE + i) = 1;
				psis[3](5 * Constants::BASIS_SIZE + i) = -1;
			}
		}
		for (auto& psi : psis)
		{
			psi.normalize();
		}

		std::vector<Resolvent_L> resolvents { 3 * NUMBER_OF_GREENSFUNCTIONS };

#pragma omp parallel for
		for (int i = 0; i < NUMBER_OF_GREENSFUNCTIONS; i++)
		{
			VectorCL a = N * psis[i];
			VectorCL b = n_hacek * psis[i];

			resolvents[3 * i].setStartingState(a);
			resolvents[3 * i + 1].setStartingState(0.5 * (a + b));
			resolvents[3 * i + 2].setStartingState(0.5 * (a + I * b));
		}
#pragma omp parallel for
		for (int i = 0; i < 3 * NUMBER_OF_GREENSFUNCTIONS; i++)
		{
			resolvents[i].compute(M, 2 * Constants::K_DISCRETIZATION);
		}

		end = std::chrono::steady_clock::now();
		std::cout << "Time for resolventes: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		return resolvents;
	}
}