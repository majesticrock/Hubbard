#include "GeneralBasis.hpp"
#include <chrono>
#include <algorithm>
#include "../MomentumIndexUtility.hpp"

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

	void GeneralBasis::printM(int i, int j) const
	{
		if (std::abs(M(j, i)) < 1e-12) {
			std::cout << 0 << "\t";
		}
		else {
			std::cout << M(j, i).real() << "\t";
		}
	}

	void GeneralBasis::printMomentumBlocks() const
	{
		std::cout << std::fixed << std::setprecision(4) << std::endl;

		for (int l = 0; l < Constants::K_DISCRETIZATION; l++)
		{
			for (int k = 0; k < Constants::K_DISCRETIZATION; ++k)
			{
				int idx = l * Constants::K_DISCRETIZATION + k;
				int jdx = addQTo(idx);
				std::cout << idx << ": " << gammaFromIndex(idx) << "\n";

				for (size_t i = 0U; i < number_of_basis_terms; ++i)
				{
					for (size_t j = 0U; j < number_of_basis_terms; ++j)
					{
						printM(j + idx * number_of_basis_terms, i + idx * number_of_basis_terms);
					}
					for (size_t j = 0U; j < number_of_basis_terms; ++j)
					{
						printM(j + jdx * number_of_basis_terms, i + idx * number_of_basis_terms);
					}
					std::cout << std::endl;
				}
				for (size_t i = 0U; i < number_of_basis_terms; ++i)
				{
					for (size_t j = 0U; j < number_of_basis_terms; ++j)
					{
						printM(j + idx * number_of_basis_terms, i + jdx * number_of_basis_terms);
					}
					for (size_t j = 0U; j < number_of_basis_terms; ++j)
					{
						printM(j + jdx * number_of_basis_terms, i + jdx * number_of_basis_terms);
					}
					std::cout << std::endl;
				}

				std::cout << std::endl;
			}
		}
	}

	void GeneralBasis::printDOSBlocks() const
	{
		std::cout << std::fixed << std::setprecision(4) << std::endl;

		for (size_t k = 0U; k < Constants::BASIS_SIZE / 2; ++k)
		{
			for (size_t i = 0U; i < number_of_basis_terms; ++i)
			{
				for (size_t j = 0U; j < number_of_basis_terms; ++j)
				{
					printM(j + k * number_of_basis_terms, i + k * number_of_basis_terms);
				}
				for (size_t j = 0U; j < number_of_basis_terms; ++j)
				{
					printM(j + (k + Constants::BASIS_SIZE / 2) * number_of_basis_terms, i + k * number_of_basis_terms);
				}
				std::cout << std::endl;
			}
			for (size_t i = 0U; i < number_of_basis_terms; ++i)
			{
				for (size_t j = 0U; j < number_of_basis_terms; ++j)
				{
					printM(j + k * number_of_basis_terms, i + (k + Constants::BASIS_SIZE / 2) * number_of_basis_terms);
				}
				for (size_t j = 0U; j < number_of_basis_terms; ++j)
				{
					printM(j + (k + Constants::BASIS_SIZE / 2) * number_of_basis_terms, i + (k + Constants::BASIS_SIZE / 2) * number_of_basis_terms);
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}

	std::vector<ResolventReturnData> GeneralBasis::computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) {
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		std::chrono::time_point end = std::chrono::steady_clock::now();

		fillMatrices();

		std::cout << "||M-M^+|| = " << (M - M.adjoint()).norm() << std::endl;
		if ((M - M.adjoint()).norm() > 1e-8) {
			throw std::runtime_error("M is not Hermitian!");
		}
		if ((N - N.adjoint()).norm() > 1e-8) {
			throw std::runtime_error("N is not Hermitian!");
		}

		end = std::chrono::steady_clock::now();
		std::cout << "Time for filling of M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		Eigen::SelfAdjointEigenSolver<MatrixCL> M_solver(M);
		Vector_L& evs_M = const_cast<Vector_L&>(M_solver.eigenvalues());
		applyMatrixOperation<OPERATION_NONE>(evs_M);
		//std::cout << "dim(kern(M)) = " << std::count_if(evs_M.begin(), evs_M.end(), [](const global_floating_type& value) { return abs(value) < 1e-16; }) << std::endl;
		//std::cout << "dim(kern(N)) = " << N.fullPivLu().kernel().cols() << std::endl;

		auto bufferMatrix = N * M_solver.eigenvectors();
		// = N * 1/M * N
		MatrixCL n_hacek = bufferMatrix
			* evs_M.unaryExpr([](global_floating_type x) { return abs(x < SALT) ? 0 : 1. / x; }).asDiagonal()
			* bufferMatrix.adjoint();

		Eigen::SelfAdjointEigenSolver<MatrixCL> norm_solver(n_hacek);
		Vector_L& evs_norm = const_cast<Vector_L&>(norm_solver.eigenvalues());
		applyMatrixOperation<OPERATION_SQRT>(evs_norm);

		// n_hacek -> n_hacek^(-1/2)
		n_hacek = norm_solver.eigenvectors()
			* evs_norm.unaryExpr([](global_floating_type x) { return abs(x < SALT) ? 0 : 1. / x; }).asDiagonal()
			* norm_solver.eigenvectors().adjoint();
		// Starting here M is the adjusted solver matrix (s s hackem)
		// n_hacek * M * n_hacek
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
		const global_floating_type norm_factor = 1. / sqrt((global_floating_type)Constants::BASIS_SIZE);
		std::vector<VectorCL> psis(NUMBER_OF_GREENSFUNCTIONS, VectorCL::Zero(TOTAL_BASIS));
		for (int i = 0; i < Constants::BASIS_SIZE; i++)
		{
			psis[0](i* number_of_basis_terms) = norm_factor;
			psis[0](i* number_of_basis_terms + 1) = norm_factor;

			psis[1](i* number_of_basis_terms) = norm_factor;
			psis[1](i* number_of_basis_terms + 1) = -norm_factor;

			if (number_of_basis_terms >= 6) {
				psis[2](i* number_of_basis_terms + 4) = norm_factor;
				psis[2](i* number_of_basis_terms + 5) = norm_factor;

				psis[3](i* number_of_basis_terms + 4) = norm_factor;
				psis[3](i* number_of_basis_terms + 5) = -norm_factor;
			}
		}

		std::vector<ResolventComplex> resolvents{ 3 * NUMBER_OF_GREENSFUNCTIONS };

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
			resolvents[i].compute(M, this->usingDOS ? 50 : 2 * Constants::K_DISCRETIZATION);
		}

		end = std::chrono::steady_clock::now();
		std::cout << "Time for resolventes: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		std::vector<ResolventReturnData> ret;
		ret.reserve(resolvents.size());
		for (const auto& re : resolvents)
		{
			ret.push_back(re.getData());
		}
		return ret;
	}
}