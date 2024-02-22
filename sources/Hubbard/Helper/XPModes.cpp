#include "XPModes.hpp"
#include <chrono>
#include <omp.h>
#include "../../Utility/PivotToBlockStructure.hpp"

namespace Hubbard::Helper {
	struct matrix_wrapper {
		Matrix_L eigenvectors;
		Vector_L eigenvalues;

		explicit matrix_wrapper(Eigen::Index size)
			: eigenvectors(Matrix_L::Zero(size, size)), eigenvalues(Vector_L::Zero(size))
		{};

		inline Matrix_L reconstruct_matrix() const
		{
			return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.adjoint();
		};
	};

	void XPModes::fillMatrices()
	{
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		constexpr int a = hermitian_size - 1;
		constexpr int b = antihermitian_size - 1;

		K_plus = Matrix_L::Zero(a * Constants::BASIS_SIZE, a * Constants::BASIS_SIZE);
		K_minus = Matrix_L::Zero(b * Constants::BASIS_SIZE, b * Constants::BASIS_SIZE);
		L = Matrix_L::Zero(a * Constants::BASIS_SIZE, b * Constants::BASIS_SIZE);

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

		const double norm_constant = this->usingDOS
			? sqrt((2.0 * this->dos_dimension) / Constants::BASIS_SIZE)
			: sqrt(1. / ((double)Constants::BASIS_SIZE));
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

		K_plus = Matrix_L::Zero(a * Constants::BASIS_SIZE, a * Constants::BASIS_SIZE);
		K_minus = Matrix_L::Zero(b * Constants::BASIS_SIZE, b * Constants::BASIS_SIZE);
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
		if ((K_plus - K_plus.adjoint()).norm() > ERROR_MARGIN * K_plus.rows() * K_plus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_plus - K_plus.adjoint()).norm()));
		if ((K_minus - K_minus.adjoint()).norm() > ERROR_MARGIN * K_minus.rows() * K_minus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_minus - K_minus.adjoint()).norm()));

		K_plus = removeNoise(K_plus);
		K_minus = removeNoise(K_minus);

		try {
			Eigen::SelfAdjointEigenSolver<Matrix_L> solver_minus(K_minus, Eigen::EigenvaluesOnly);
			Vector_L& evs = const_cast<Vector_L&>(solver_minus.eigenvalues());
			applyMatrixOperation<OPERATION_NONE>(evs);

			Eigen::SelfAdjointEigenSolver<Matrix_L> solver_plus(K_plus, Eigen::EigenvaluesOnly);
			evs = const_cast<Vector_L&>(solver_plus.eigenvalues());
			applyMatrixOperation<OPERATION_NONE>(evs);
		}
		catch (const MatrixIsNegativeException& ex) {
			return true;
		}

		return false;
	};

	std::vector<ResolventReturnData> XPModes::computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever)
	{
		fillMatrices();
		createStartingStates();

		// M_new = K_plus
		Matrix_L solver_matrix;
		matrix_wrapper k_solutions[2] = { matrix_wrapper(K_plus.rows()), matrix_wrapper(K_minus.rows()) };

		omp_set_nested(2);
		Eigen::initParallel();

#pragma omp parallel sections
		{
#pragma omp section
			{
				std::chrono::time_point begin_in = std::chrono::steady_clock::now();
				auto pivot = Utility::pivot_to_block_structure(K_plus);
				K_plus = pivot.transpose() * K_plus * pivot;
				auto blocks = Utility::identify_hermitian_blocks(K_plus);

				std::chrono::time_point end_in = std::chrono::steady_clock::now();
				std::cout << "Time for pivoting K_+: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
				begin_in = std::chrono::steady_clock::now();

#pragma omp parallel for
				for (int i = 0; i < blocks.size(); ++i)
				{
					Eigen::SelfAdjointEigenSolver<Matrix_L> solver(K_plus.block(blocks[i].first, blocks[i].first, blocks[i].second, blocks[i].second));
					k_solutions[0].eigenvalues.segment(blocks[i].first, blocks[i].second) = solver.eigenvalues();
					k_solutions[0].eigenvectors.block(blocks[i].first, blocks[i].first, blocks[i].second, blocks[i].second) = solver.eigenvectors();
				}

				applyMatrixOperation<OPERATION_NONE>(k_solutions[0].eigenvalues);
				k_solutions[0].eigenvectors.applyOnTheLeft(pivot);

				end_in = std::chrono::steady_clock::now();
				std::cout << "Time for solving K_+: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
			}
#pragma omp section
			{
				std::chrono::time_point begin_in = std::chrono::steady_clock::now();

				auto pivot = Utility::pivot_to_block_structure(K_minus);
				K_minus = pivot.transpose() * K_minus * pivot;
				auto blocks = Utility::identify_hermitian_blocks(K_minus);

				std::chrono::time_point end_in = std::chrono::steady_clock::now();
				std::cout << "Time for pivoting K_-: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
				begin_in = std::chrono::steady_clock::now();

#pragma omp parallel for
				for (int i = 0; i < blocks.size(); ++i)
				{
					Eigen::SelfAdjointEigenSolver<Matrix_L> solver(K_minus.block(blocks[i].first, blocks[i].first, blocks[i].second, blocks[i].second));
					k_solutions[1].eigenvalues.segment(blocks[i].first, blocks[i].second) = solver.eigenvalues();
					k_solutions[1].eigenvectors.block(blocks[i].first, blocks[i].first, blocks[i].second, blocks[i].second) = solver.eigenvectors();
				}

				applyMatrixOperation<OPERATION_NONE>(k_solutions[1].eigenvalues);
				k_solutions[1].eigenvectors.applyOnTheLeft(pivot);

				end_in = std::chrono::steady_clock::now();
				std::cout << "Time for solving K_-: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
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

			//Eigen::SelfAdjointEigenSolver<Matrix_L> solver;
			//solver.compute(N_new);

			auto pivot = Utility::pivot_to_block_structure(N_new);
			N_new = pivot.transpose() * N_new * pivot;
			auto blocks = Utility::identify_hermitian_blocks(N_new);

			matrix_wrapper n_solution{ N_new.rows() };
#pragma omp parallel for
			for (int i = 0; i < blocks.size(); ++i)
			{
				Eigen::SelfAdjointEigenSolver<Matrix_L> solver(N_new.block(blocks[i].first, blocks[i].first, blocks[i].second, blocks[i].second));
				n_solution.eigenvalues.segment(blocks[i].first, blocks[i].second) = solver.eigenvalues();
				n_solution.eigenvectors.block(blocks[i].first, blocks[i].first, blocks[i].second, blocks[i].second) = solver.eigenvectors();
			}
			n_solution.eigenvectors.applyOnTheLeft(pivot);

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