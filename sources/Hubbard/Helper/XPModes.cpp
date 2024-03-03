#include "XPModes.hpp"
#include <chrono>
#include <omp.h>
#include "../../Utility/PivotToBlockStructure.hpp"
#include <Eigen/Sparse>

namespace Hubbard::Helper {
	XPModes::matrix_wrapper XPModes::matrix_wrapper::pivot_and_solve(Matrix_L& toSolve)
	{
		auto pivot = Utility::pivot_to_block_structure(toSolve);
		toSolve = pivot.transpose() * toSolve * pivot;
		auto blocks = Utility::identify_hermitian_blocks(toSolve);
		matrix_wrapper solution(toSolve.rows());

#pragma omp parallel for
		for (int i = 0; i < blocks.size(); ++i)
		{
			Eigen::SelfAdjointEigenSolver<Matrix_L> solver(toSolve.block(blocks[i].position, blocks[i].position, blocks[i].size, blocks[i].size));
			solution.eigenvalues.segment(blocks[i].position, blocks[i].size) = solver.eigenvalues();
			solution.eigenvectors.block(blocks[i].position, blocks[i].position, blocks[i].size, blocks[i].size) = solver.eigenvectors();
		}
		solution.eigenvectors.applyOnTheLeft(pivot);
		return solution;
	}
	bool XPModes::matrix_wrapper::is_non_negative(Matrix_L& toSolve)
	{
		auto pivot = Utility::pivot_to_block_structure(toSolve);
		toSolve = pivot.transpose() * toSolve * pivot;
		auto blocks = Utility::identify_hermitian_blocks(toSolve);
		for (const auto& block : blocks)
		{
			Eigen::SelfAdjointEigenSolver<Matrix_L> solver(toSolve.block(block.position, block.position, block.size, block.size), Eigen::EigenvaluesOnly);
			//Eigen::LDLT<Matrix_L> cholesky(toSolve.block(block.position, block.position, block.size, block.size));
			if (ModeHelper::contains_negative(solver.eigenvalues())) {
				return false;
			}
		}
		return true;
	};

	void XPModes::fillMatrices()
	{
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		constexpr int a = hermitian_size - 1;
		constexpr int b = antihermitian_size - 1;

		K_plus.setZero(a * Constants::BASIS_SIZE, a * Constants::BASIS_SIZE);
		K_minus.setZero(b * Constants::BASIS_SIZE, b * Constants::BASIS_SIZE);
		L.setZero(a * Constants::BASIS_SIZE, b * Constants::BASIS_SIZE);

		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size)) {
					fill_block_M(i, j);
				}
				// N only contains offdiagonal blocks
				else if (i < hermitian_size && j >= hermitian_size) {
					fill_block_N(i, j);
				}
			}
		}

		if ((K_plus - K_plus.adjoint()).norm() > ERROR_MARGIN * K_plus.rows() * K_plus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_plus - K_plus.adjoint()).norm()));
		if ((K_minus - K_minus.adjoint()).norm() > ERROR_MARGIN * K_minus.rows() * K_minus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_minus - K_minus.adjoint()).norm()));

		L = removeNoise(L);
		K_plus = removeNoise(K_plus);
		K_minus = removeNoise(K_minus);

		std::chrono::time_point end = std::chrono::steady_clock::now();
		std::cout << "Time for filling of M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	}

	void XPModes::createStartingStates()
	{
		this->startingState_SC = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };
		this->startingState_CDW = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };
		this->startingState_AFM = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };
		this->startingState_AFM_transversal = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };

		const global_floating_type norm_constant =
#ifdef _EXACT_DOS
			sqrt(1. / ((global_floating_type)Constants::BASIS_SIZE));
#else
			this->usingDOS ? sqrt((2.0 * this->dos_dimension) / Constants::BASIS_SIZE) : sqrt(1. / ((global_floating_type)Constants::BASIS_SIZE));
#endif
		for (int j = 0; j < 2; ++j)
		{
			for (size_t i = 0U; i < Constants::BASIS_SIZE; ++i)
			{
				startingState_SC[j](i) = norm_constant;
				startingState_CDW[j](2 * Constants::BASIS_SIZE + i) = norm_constant;
				startingState_AFM[j](2 * Constants::BASIS_SIZE + i) = (i < Constants::BASIS_SIZE / 2) ? norm_constant : -norm_constant;
				startingState_AFM_transversal[j](3 * Constants::BASIS_SIZE + i) = norm_constant;
			}
		}
	}

	bool XPModes::matrix_is_negative() {
		constexpr int a = hermitian_size - 1;
		constexpr int b = antihermitian_size - 1;

		K_plus.setZero(a * Constants::BASIS_SIZE, a * Constants::BASIS_SIZE);
		K_minus.setZero(b * Constants::BASIS_SIZE, b * Constants::BASIS_SIZE);
		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size)) {
					fill_block_M(i, j);
				}
			}
		}
#ifdef _DEBUG
		if ((K_plus - K_plus.adjoint()).norm() > ERROR_MARGIN * K_plus.rows() * K_plus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_plus - K_plus.adjoint()).norm()));
		if ((K_minus - K_minus.adjoint()).norm() > ERROR_MARGIN * K_minus.rows() * K_minus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_minus - K_minus.adjoint()).norm()));
#endif
		if (contains_negative(K_minus.diagonal()) || contains_negative(K_plus.diagonal())) {
			return true;
		}
		K_minus = removeNoise(K_minus);
		if (not matrix_wrapper::is_non_negative(K_minus)) {
			return true;
		}
		K_plus = removeNoise(K_plus);
		if (not matrix_wrapper::is_non_negative(K_plus)) {
			return true;
		}
		return false;
	};

	std::vector<ResolventReturnData> XPModes::computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever)
	{
		fillMatrices();
		createStartingStates();

		Matrix_L solver_matrix;
		matrix_wrapper k_solutions[2];

		omp_set_nested(2);
		Eigen::initParallel();

#pragma omp parallel sections
		{
#pragma omp section
			{
				std::chrono::time_point begin_in = std::chrono::steady_clock::now();
				k_solutions[0] = matrix_wrapper::pivot_and_solve(K_plus);
				applyMatrixOperation<OPERATION_NONE>(k_solutions[0].eigenvalues);
				std::chrono::time_point end_in = std::chrono::steady_clock::now();
				std::cout << "Time for solving K_+: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;

				// free the allocated memory
				K_plus.conservativeResize(0, 0);
			}
#pragma omp section
			{
				std::chrono::time_point begin_in = std::chrono::steady_clock::now();
				k_solutions[1] = matrix_wrapper::pivot_and_solve(K_minus);
				applyMatrixOperation<OPERATION_NONE>(k_solutions[1].eigenvalues);
				std::chrono::time_point end_in = std::chrono::steady_clock::now();
				std::cout << "Time for solving K_-: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;

				// free the allocated memory
				K_minus.conservativeResize(0, 0);
			}
		}

		/* plus(minus)_index indicates whether the upper left block is for the
		* Hermitian or the anti-Hermitian operators.
		* The default is that the upper left block contains the Hermtian operators,
		* then plus_index = 0 and minus_index = 1
		*/
		auto compute_solver_matrix = [&](size_t plus_index, size_t minus_index) {
			std::chrono::time_point begin_in = std::chrono::steady_clock::now();
			if (minus_index == 0) L.transposeInPlace();
			solver_matrix.resize(k_solutions[plus_index].eigenvalues.rows(), k_solutions[plus_index].eigenvalues.rows());

			Vector_L K_EV = k_solutions[minus_index].eigenvalues;
			applyMatrixOperation<OPERATION_INVERSE>(K_EV);
			Matrix_L buffer_matrix = L * k_solutions[minus_index].eigenvectors;
			Matrix_L N_new = buffer_matrix * K_EV.asDiagonal() * buffer_matrix.adjoint();

			std::chrono::time_point end_in = std::chrono::steady_clock::now();
			std::cout << "Time for computing N_new: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
			begin_in = std::chrono::steady_clock::now();

			matrix_wrapper n_solution = matrix_wrapper::pivot_and_solve(N_new);
			applyMatrixOperation<OPERATION_INVERSE_SQRT>(n_solution.eigenvalues);

			// Starting here, N_new = 1/sqrt(N_new)
			// I forego another matrix to save some memory
			N_new = n_solution.eigenvectors * n_solution.eigenvalues.asDiagonal() * n_solution.eigenvectors.adjoint();
			startingState_SC[plus_index] = N_new * L * startingState_SC[plus_index];
			startingState_CDW[plus_index] = N_new * L * startingState_CDW[plus_index];
			startingState_AFM[plus_index] = N_new * L * startingState_AFM[plus_index];
			startingState_AFM_transversal[plus_index] = N_new * L * startingState_AFM_transversal[plus_index];

			end_in = std::chrono::steady_clock::now();
			std::cout << "Time for adjusting N_new: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;

			begin_in = std::chrono::steady_clock::now();
			buffer_matrix = N_new * k_solutions[plus_index].eigenvectors;
			solver_matrix = removeNoise((buffer_matrix * k_solutions[plus_index].eigenvalues.asDiagonal() * buffer_matrix.adjoint()).eval());
			end_in = std::chrono::steady_clock::now();
			std::cout << "Time for computing solver_matrix: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
			}; // end lambda

		std::chrono::time_point begin = std::chrono::steady_clock::now();

		constexpr unsigned int N_RESOLVENT_TYPES = 4U;
		std::vector<ResolventReal> resolvents{};
		resolvents.reserve(2 * N_RESOLVENT_TYPES);

		int LANCZOS_ITERATION_NUMBER;
		if (this->usingDOS) {
			LANCZOS_ITERATION_NUMBER = 150;
		}
		else {
			LANCZOS_ITERATION_NUMBER = 2 * Constants::K_DISCRETIZATION;
		}

		for (size_t i = 0U; i < 2U; ++i)
		{
			// It is going to compute the anti-Hermitian first
			compute_solver_matrix(i, 1 - i);
			resolvents.push_back(ResolventReal(startingState_SC[i]));
			resolvents.push_back(ResolventReal(startingState_CDW[i]));
			resolvents.push_back(ResolventReal(startingState_AFM[i]));
			resolvents.push_back(ResolventReal(startingState_AFM_transversal[i]));

#pragma omp parallel sections
			{
#pragma omp section
				{
					resolvents[N_RESOLVENT_TYPES * i].compute(solver_matrix, LANCZOS_ITERATION_NUMBER);
				}
#pragma omp section
				{
					resolvents[N_RESOLVENT_TYPES * i + 1].compute(solver_matrix, LANCZOS_ITERATION_NUMBER);
				}
#pragma omp section
				{
					resolvents[N_RESOLVENT_TYPES * i + 2].compute(solver_matrix, LANCZOS_ITERATION_NUMBER);
				}
#pragma omp section
				{
					resolvents[N_RESOLVENT_TYPES * i + 3].compute(solver_matrix, LANCZOS_ITERATION_NUMBER);
				}
			}
		}

		std::chrono::time_point end = std::chrono::steady_clock::now();
		std::cout << "Time for resolvents: "
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