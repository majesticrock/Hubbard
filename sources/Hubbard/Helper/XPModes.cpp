#include "XPModes.hpp"
#include <chrono>
#include <omp.h>

namespace Hubbard::Helper {
	void XPModes::fillMatrices()
	{
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

		auto setZero = [](global_floating_type val) {
			return (abs(val) < DEFAULT_PRECISION ? 0 : val);
			};
		K_plus = K_plus.array().unaryExpr(setZero);
		K_minus = K_minus.array().unaryExpr(setZero);

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
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		std::chrono::time_point end = std::chrono::steady_clock::now();
		
		fillMatrices();

		if ((K_plus - K_plus.adjoint()).norm() > ERROR_MARGIN * K_plus.rows() * K_plus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_plus - K_plus.adjoint()).norm()));
		if ((K_minus - K_minus.adjoint()).norm() > ERROR_MARGIN * K_minus.rows() * K_minus.cols())
			throw std::invalid_argument("K_+ is not hermitian: " + to_string((K_minus - K_minus.adjoint()).norm()));

		auto setZero = [](global_floating_type val) {
			return (abs(val) < DEFAULT_PRECISION ? 0 : val);
			};
		L = L.array().unaryExpr(setZero);
		K_plus = K_plus.array().unaryExpr(setZero);
		K_minus = K_minus.array().unaryExpr(setZero);

		end = std::chrono::steady_clock::now();
		std::cout << "Time for filling of M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		Vector_L startingState_SC[2] = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };
		Vector_L startingState_CDW[2] = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };
		Vector_L startingState_AFM[2] = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };
		Vector_L startingState_AFM_transversal[2] = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };

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

		// M_new = K_plus
		Matrix_L solver_matrix;
		Eigen::SelfAdjointEigenSolver<Matrix_L> k_solver[2];

		omp_set_nested(1);
		Eigen::initParallel();

#pragma omp parallel sections
		{
#pragma omp section
			{
				std::chrono::time_point begin_in = std::chrono::steady_clock::now();
				k_solver[0].compute(K_plus);

				global_floating_type tol{ k_solver[0].eigenvalues().norm() * ERROR_MARGIN };
				if (k_solver[0].info() == Eigen::Success && (k_solver[0].eigenvalues().array() >= -tol).all()) {
					std::cout << "K_+: eigenvalues are real and accurate up to tolerance" << std::endl;
				}
				else {
					std::cerr << "K_+: eigenvalues may not be accurate: " << k_solver[0].info() << std::endl;
				}

				// Do the most scary thing I've ever seen in c++
				// And remove the const from the reference returned by eigenvalues()
				Vector_L& evs = const_cast<Vector_L&>(k_solver[0].eigenvalues());
				applyMatrixOperation<OPERATION_NONE>(evs);

				std::chrono::time_point end_in = std::chrono::steady_clock::now();
				std::cout << "Time for solving K_+: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
			}
#pragma omp section
			{
				std::chrono::time_point begin_in = std::chrono::steady_clock::now();
				k_solver[1].compute(K_minus);

				global_floating_type tol{ k_solver[1].eigenvalues().norm() * ERROR_MARGIN };
				if (k_solver[1].info() == Eigen::Success && (k_solver[1].eigenvalues().array() >= -tol).all()) {
					std::cout << "K_-: eigenvalues are real and accurate up to tolerance" << std::endl;
				}
				else {
					std::cerr << "K_-: eigenvalues may not be accurate: " << k_solver[1].info() << std::endl;
				}

				// Do the most scary thing I've ever seen in c++
				// And remove the const from the reference returned by eigenvalues()
				Vector_L& evs = const_cast<Vector_L&>(k_solver[1].eigenvalues());
				applyMatrixOperation<OPERATION_NONE>(evs);

				std::chrono::time_point end_in = std::chrono::steady_clock::now();
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
			solver_matrix.resize(k_solver[plus_index].eigenvalues().rows(), k_solver[plus_index].eigenvalues().rows());

			Vector_L K_EV = k_solver[minus_index].eigenvalues();
			applyMatrixOperation<OPERATION_INVERSE>(K_EV);
			Matrix_L buffer_matrix = L * k_solver[minus_index].eigenvectors();
			Matrix_L N_new = buffer_matrix * K_EV.asDiagonal() * buffer_matrix.adjoint();

			std::chrono::time_point end_in = std::chrono::steady_clock::now();
			std::cout << "Time for computing N_new: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
			begin_in = std::chrono::steady_clock::now();

			Eigen::SelfAdjointEigenSolver<Matrix_L> solver;
			solver.compute(N_new);
			global_floating_type tol{ solver.eigenvalues().norm() * ERROR_MARGIN };
			if (solver.info() == Eigen::Success && (solver.eigenvalues().array() >= -tol).all()) {
				std::cout << "N_new: eigenvalues are real and accurate up to tolerance" << std::endl;
			}
			else {
				std::cerr << "N_new: eigenvalues may not be accurate: " << solver.info() << std::endl;
			}

			Vector_L& n_ev = const_cast<Vector_L&>(solver.eigenvalues());
			applyMatrixOperation<OPERATION_INVERSE_SQRT>(n_ev);

			// Starting here, N_new = 1/sqrt(N_new)
			// I forego another matrix to save some memory
			N_new = solver.eigenvectors() * n_ev.asDiagonal() * solver.eigenvectors().adjoint();
			startingState_SC[plus_index] = N_new * L * startingState_SC[plus_index];
			startingState_CDW[plus_index] = N_new * L * startingState_CDW[plus_index];
			startingState_AFM[plus_index] = N_new * L * startingState_AFM[plus_index];
			startingState_AFM_transversal[plus_index] = N_new * L * startingState_AFM_transversal[plus_index];

			end_in = std::chrono::steady_clock::now();
			std::cout << "Time for adjusting N_new: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;

			begin_in = std::chrono::steady_clock::now();
			buffer_matrix = N_new * k_solver[plus_index].eigenvectors();
			solver_matrix = (buffer_matrix * k_solver[plus_index].eigenvalues().asDiagonal() * buffer_matrix.adjoint()).array().unaryExpr(setZero);
			end_in = std::chrono::steady_clock::now();
			std::cout << "Time for computing solver_matrix: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
			}; // end lambda

		begin = std::chrono::steady_clock::now();

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

		end = std::chrono::steady_clock::now();
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