#include "GeneralBasis.hpp"
#include <chrono>

namespace Hubbard::Helper {
	void GeneralBasis::fillBlock(int i, int j) {
		// fill N
		for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0) {
					if (term.delta_momenta.size() > 1) throw std::invalid_argument("Too many deltas: " + term.delta_momenta.size());
					if (term.delta_momenta[0].first.momentum_list.size() != 1) throw std::invalid_argument("First delta list is not of size 1: " + term.delta_momenta[0].first.momentum_list.size());
					if (term.delta_momenta[0].second.momentum_list.size() != 1) throw std::invalid_argument("To be implemented: Second delta list is not of size 1: " + term.delta_momenta[0].second.momentum_list.size());

					int l_buf = k;
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						Eigen::Vector2i l_buf_vec = { x(k), y(k) };
						l_buf_vec(0) += Constants::K_DISCRETIZATION;
						l_buf_vec(1) += Constants::K_DISCRETIZATION;
						clean_factor_2pi(l_buf_vec);
						l_buf = l_buf_vec(0) * 2 * Constants::K_DISCRETIZATION + l_buf_vec(1);
					}
					N(j * Constants::BASIS_SIZE + l_buf, i * Constants::BASIS_SIZE + k) += computeTerm(term, l_buf, k);
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; l++)
					{
						N(j * Constants::BASIS_SIZE + l, i * Constants::BASIS_SIZE + k) += computeTerm(term, l, k);
					}
				}
			}
		}

		// fill M
		for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0) {
					if (term.delta_momenta.size() > 1) throw std::invalid_argument("Too many deltas: " + term.delta_momenta.size());
					if (term.delta_momenta[0].first.momentum_list.size() != 1) throw std::invalid_argument("First delta list is not of size 1: " + term.delta_momenta[0].first.momentum_list.size());
					if (term.delta_momenta[0].second.momentum_list.size() != 1) throw std::invalid_argument("To be implemented: Second delta list is not of size 1: " + term.delta_momenta[0].second.momentum_list.size());

					int l_buf = k;
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						Eigen::Vector2i l_buf_vec = { x(k), y(k) };
						l_buf_vec(0) += Constants::K_DISCRETIZATION;
						l_buf_vec(1) += Constants::K_DISCRETIZATION;
						clean_factor_2pi(l_buf_vec);
						l_buf = l_buf_vec(0) * 2 * Constants::K_DISCRETIZATION + l_buf_vec(1);
					}
					M(j * Constants::BASIS_SIZE + l_buf, i * Constants::BASIS_SIZE + k) += computeTerm(term, l_buf, k);
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; l++)
					{
						M(j * Constants::BASIS_SIZE + l, i * Constants::BASIS_SIZE + k) += computeTerm(term, l, k);
					}
				}
			}
		}
	}

	void GeneralBasis::fillMatrices()
	{
		M = MatrixCL::Zero(TOTAL_BASIS, TOTAL_BASIS);
		N = MatrixCL::Zero(TOTAL_BASIS, TOTAL_BASIS);

		//#pragma omp parallel for
		for (int i = 0; i < number_of_basis_terms; i++)
		{
			for (int j = 0; j < number_of_basis_terms; j++)
			{
				fillBlock(i, j);
			}
		}
	}

	std::vector<Resolvent_L> GeneralBasis::computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) {
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		std::chrono::time_point end = std::chrono::steady_clock::now();

		fillMatrices();

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

			psis[2](4 * Constants::BASIS_SIZE + i) = 1;
			psis[2](5 * Constants::BASIS_SIZE + i) = 1;

			psis[3](4 * Constants::BASIS_SIZE + i) = 1;
			psis[3](5 * Constants::BASIS_SIZE + i) = -1;
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