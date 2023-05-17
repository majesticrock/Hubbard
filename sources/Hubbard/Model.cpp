#include "Model.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <fstream>
#include <algorithm>
#include "../Utility/Lanczos.hpp"

// Both methods yield precisely the same data!
#define _PSEUDO_INVERSE

namespace Hubbard {
	constexpr double_prec SQRT_SALT = 1e-5;
	constexpr double_prec SALT = SQRT_SALT * SQRT_SALT;
	constexpr double_prec ERROR_MARGIN = 1e-10;

	void Model::initializeParameters()
	{
		this->BASIS_SIZE = 4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION;
		if (start_basis_at < 0) {
			// We investigate the special x-p-basis
			this->TOTAL_BASIS = this->BASIS_SIZE * 8;
			this->number_of_basis_terms = 10;
		}
		else {
			this->TOTAL_BASIS = this->BASIS_SIZE * this->number_of_basis_terms;
		}
		this->delta_cdw_up = 0.1;
		this->delta_cdw_down = (U > 0) ? -this->delta_cdw_up : this->delta_cdw_up;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
		computeChemicalPotential();
	}

	Model::Model(double_prec _temperature, double_prec _U, int _number_of_basis_terms, int _start_basis_at)
		: temperature(_temperature), U(_U), number_of_basis_terms(_number_of_basis_terms), start_basis_at(_start_basis_at)
	{
		initializeParameters();
	}

	Model::Model(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: temperature(_params.temperature), U(_params.U), number_of_basis_terms(_number_of_basis_terms), start_basis_at(_start_basis_at)
	{
		initializeParameters();
	}

	std::unique_ptr<std::vector<Resolvent_L>> Model::computeCollectiveModes(std::vector<std::vector<double>>& reciever)
	{
		std::cout << "Gap values:  " << delta_cdw_up << "  " << delta_cdw_down << "  " << delta_sc << "  "
			<< delta_eta << " " << delta_occupation_up << " " << delta_occupation_down << std::endl;
		// First off we need to compute every possible expectation value
		// We use the mean field system's symmetries
		// i.e. there are only the standard SC, CDW, Eta and N operators non-zero
		// spin symmetry and conservation of momentum up to addition of Q holds
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		std::chrono::time_point end = std::chrono::steady_clock::now();

		{
			Matrix_L rho = Matrix_L::Zero(4, 4);
			Eigen::SelfAdjointEigenSolver<Matrix_L> solver;

			expecs = std::vector<Matrix_L>(4, Matrix_L::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
			for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
			{
				for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
				{
					fillHamiltonian((k * L_PI) / Constants::K_DISCRETIZATION, (l * L_PI) / Constants::K_DISCRETIZATION);
					solver.compute(hamilton);
					rho.fill(0);
					for (int i = 0; i < 4; i++)
					{
						rho(i, i) = fermi_dirac(solver.eigenvalues()(i));
					}
					rho = solver.eigenvectors() * rho * (solver.eigenvectors().transpose());
					for (int idx = 0; idx < 4; idx++)
					{
						expecs[idx](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION) = rho(idx, 0);
						sum_of_all[idx] += rho(idx, 0);
					}
					if (std::abs(rho(3, 0)) > SALT) {
						std::cerr << "Warning: <eta> does not vanish! " << rho(3, 0) << std::endl;
					}
				}
			}
		}
		computeChemicalPotential();
		loadWick("../commutators/wick_");
		end = std::chrono::steady_clock::now();
		std::cout << "Time for expectation values: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		if (this->start_basis_at < 0) {
			fill_M_N_xp_basis();
		}
		else {
			fill_M_N();
			Eigen::SelfAdjointEigenSolver<Matrix_L> solver(M);
			for (size_t i = 0; i < solver.eigenvalues().size(); i++)
			{
				if (solver.eigenvalues()(i) < 0) {
					std::cout << solver.eigenvalues()(i) << "\n\n" << solver.eigenvectors().col(i) << std::endl << std::endl;
				}
			}

			std::cout << M << std::endl;
			return nullptr;
		}

		{
			for (size_t i = 0; i < L.rows(); i++)
			{
				for (size_t j = 0; j < L.cols(); j++)
				{
					if (std::abs(L(i, j)) < SALT) {
						L(i, j) = 0;
					}
				}
			}
			for (size_t i = 0; i < K_plus.rows(); i++)
			{
				for (size_t j = 0; j < K_plus.cols(); j++)
				{
					if (std::abs(K_plus(i, j)) < SALT) {
						K_plus(i, j) = 0;
					}
					if (std::abs(K_plus(i, j) - K_plus(j, i)) > ERROR_MARGIN) {
						std::cerr << std::scientific << std::setprecision(12) << std::abs(K_plus(i, j) - K_plus(j, i)) << std::endl;
						throw std::invalid_argument("K_+ is not hermitian: " + std::to_string(K_plus(i, j))
							+ " || " + std::to_string(K_plus(j, i)) + "\t\tPosition: " + std::to_string(i) + ", " + std::to_string(j));
					}
				}
			}
			for (size_t i = 0; i < K_minus.rows(); i++)
			{
				for (size_t j = 0; j < K_minus.cols(); j++)
				{
					if (std::abs(K_minus(i, j)) < SALT) {
						K_minus(i, j) = 0;
					}
					if (std::abs(K_minus(i, j) - K_minus(j, i)) > ERROR_MARGIN) {
						std::cerr << std::scientific << std::setprecision(12) << std::abs(K_minus(i, j) - K_minus(j, i)) << std::endl;
						throw std::invalid_argument("K_- is not hermitian: " + std::to_string(K_minus(i, j))
							+ " || " + std::to_string(K_minus(j, i)) + "\t\tPosition: " + std::to_string(i) + ", " + std::to_string(j));
					}
				}
			}
		}

		end = std::chrono::steady_clock::now();
		std::cout << "Time for filling of M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		Vector_L startingState_SC[2] = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };
		Vector_L startingState_CDW[2] = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };
		Vector_L startingState_AFM[2] = { Vector_L::Zero(K_minus.rows()),  Vector_L::Zero(K_plus.rows()) };
		for (int j = 0; j < 2; j++)
		{
			for (size_t i = 0; i < BASIS_SIZE; i++)
			{
				startingState_SC[j](i) = 1;
				startingState_CDW[j](2 * BASIS_SIZE + i) = 1;
				startingState_AFM[j](2 * BASIS_SIZE + i) = (i < BASIS_SIZE / 2) ? 1 : -1;
			}
			startingState_SC[j].normalize();
			startingState_CDW[j].normalize();
			startingState_AFM[j].normalize();
		}

		// M_new = K_plus
		Matrix_L solver_matrix;
		Eigen::SelfAdjointEigenSolver<Matrix_L> k_solver[2];

		omp_set_num_threads(8);
		omp_set_nested(1);
		Eigen::initParallel();

#pragma omp parallel sections
		{
#pragma omp section
			{
				std::chrono::time_point begin_in = std::chrono::steady_clock::now();
				k_solver[0].compute(K_plus);

				double_prec tol = k_solver[0].eigenvalues().norm() * ERROR_MARGIN;
				if (k_solver[0].info() == Eigen::Success && (k_solver[0].eigenvalues().array() >= -tol).all()) {
					std::cout << "K_+: eigenvalues are real and accurate up to tolerance" << std::endl;
				}
				else {
					std::cerr << "K_+: eigenvalues may not be accurate: " << k_solver[0].info() << std::endl;
				}

				// Do the most scary thing I've ever seen in c++
				// And remove the const from the reference returned by eigenvalues()
				Vector_L& evs = const_cast<Vector_L&>(k_solver[0].eigenvalues());
				for (size_t i = 0; i < evs.size(); i++)
				{
					if (evs(i) < -SALT) {
						std::cerr << "K_+:   " << evs(i) << std::endl;
						//throw std::invalid_argument("K_+ is not positive!  " + std::to_string(evs(i)));
					}
					if (evs(i) < SALT) {
#ifdef _PSEUDO_INVERSE
						evs(i) = 0;
#else
						evs(i) = SALT;
#endif
					}
				}

				std::chrono::time_point end_in = std::chrono::steady_clock::now();
				std::cout << "Time for solving K_+: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
			}
#pragma omp section
			{
				std::chrono::time_point begin_in = std::chrono::steady_clock::now();
				k_solver[1].compute(K_minus);

				std::cout << "K_- recon:    "
					<< (k_solver[1].eigenvectors() * k_solver[1].eigenvalues().asDiagonal() * k_solver[1].eigenvectors().transpose()
						- K_minus).norm() << std::endl;

				double_prec tol = k_solver[1].eigenvalues().norm() * ERROR_MARGIN;
				if (k_solver[1].info() == Eigen::Success && (k_solver[1].eigenvalues().array() >= -tol).all()) {
					std::cout << "K_-: eigenvalues are real and accurate up to tolerance" << std::endl;
				}
				else {
					std::cerr << "K_-: eigenvalues may not be accurate: " << k_solver[1].info() << std::endl;
				}

				// Do the most scary thing I've ever seen in c++
				// And remove the const from the reference returned by eigenvalues()
				Vector_L& evs = const_cast<Vector_L&>(k_solver[1].eigenvalues());
				for (size_t i = 0; i < evs.size(); i++)
				{
					if (evs(i) < -SALT) {
						std::cerr << "K_-:   " << evs(i) << std::endl;
						//throw std::invalid_argument("K_- is not positive!  " + std::to_string(evs(i)));
					}
					if (evs(i) < SALT) {
#ifdef _PSEUDO_INVERSE
						evs(i) = 0;
#else
						evs(i) = SALT;
#endif
					}
				}

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
			for (size_t i = 0; i < k_solver[minus_index].eigenvalues().size(); i++)
			{
				if (K_EV(i) > SALT) {
					K_EV(i) = 1. / (k_solver[minus_index].eigenvalues()(i));
				}
				else {
#ifdef _PSEUDO_INVERSE
					K_EV(i) = 0;
#else
					K_EV(i) = SALT;
#endif
				}
			}
			Matrix_L buffer_matrix = L * k_solver[minus_index].eigenvectors();
			Matrix_L N_new = buffer_matrix * K_EV.asDiagonal() * buffer_matrix.transpose();

			std::chrono::time_point end_in = std::chrono::steady_clock::now();
			std::cout << "Time for computing N_new: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
			begin_in = std::chrono::steady_clock::now();

			Eigen::SelfAdjointEigenSolver<Matrix_L> solver;
			solver.compute(N_new);
			double_prec tol = solver.eigenvalues().norm() * ERROR_MARGIN;
			if (solver.info() == Eigen::Success && (solver.eigenvalues().array() >= -tol).all()) {
				std::cout << "N_new: eigenvalues are real and accurate up to tolerance" << std::endl;
			}
			else {
				std::cerr << "N_new: eigenvalues may not be accurate: " << solver.info() << std::endl;
			}

			Vector_L& n_ev = const_cast<Vector_L&>(solver.eigenvalues());
			for (size_t i = 0; i < n_ev.size(); i++)
			{
				if (solver.eigenvalues()(i) < -ERROR_MARGIN) {
					std::cerr << solver.eigenvalues()(i) << std::endl;
					throw std::invalid_argument("N_new is not positive!  " + std::to_string(solver.eigenvalues()(i)));
				}
				else if (n_ev(i) < SALT) {
#ifdef _PSEUDO_INVERSE
					n_ev(i) = 0;
#else
					n_ev(i) = SQRT_SALT;
#endif
				}
				else {
					n_ev(i) = 1. / sqrt(n_ev(i));
				}
			}
			// Starting here, N_new = 1/sqrt(N_new)
			// I forego another matrix to save some memory
			N_new = solver.eigenvectors() * n_ev.asDiagonal() * solver.eigenvectors().adjoint();
			startingState_SC[plus_index] = N_new * L * startingState_SC[plus_index];
			startingState_CDW[plus_index] = N_new * L * startingState_CDW[plus_index];
			startingState_AFM[plus_index] = N_new * L * startingState_AFM[plus_index];

			end_in = std::chrono::steady_clock::now();
			std::cout << "Time for adjusting N_new: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;

			begin_in = std::chrono::steady_clock::now();
			buffer_matrix = N_new * k_solver[plus_index].eigenvectors();
			solver_matrix = (buffer_matrix * k_solver[plus_index].eigenvalues().asDiagonal() * buffer_matrix.adjoint());
			end_in = std::chrono::steady_clock::now();
			std::cout << "Time for computing solver_matrix: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
		}; // end lambda

		begin = std::chrono::steady_clock::now();

		std::unique_ptr<std::vector<Resolvent_L>> resolvents = std::make_unique< std::vector<Resolvent_L>>(std::vector<Resolvent_L>());
		resolvents->reserve(6);

		for (size_t i = 0; i < 2; i++)
		{
			// It is going to compute the anti-Hermitian first
			compute_solver_matrix(i, 1 - i);
			resolvents->push_back(Resolvent_L(startingState_SC[i]));
			resolvents->push_back(Resolvent_L(startingState_CDW[i]));
			resolvents->push_back(Resolvent_L(startingState_AFM[i]));
#pragma omp parallel sections
			{
#pragma omp section
				{
					(*resolvents)[3 * i].compute(solver_matrix, 2 * Constants::K_DISCRETIZATION);
				}
#pragma omp section
				{
					(*resolvents)[3 * i + 1].compute(solver_matrix, 2 * Constants::K_DISCRETIZATION);
				}
#pragma omp section
				{
					(*resolvents)[3 * i + 2].compute(solver_matrix, 2 * Constants::K_DISCRETIZATION);
				}
			}
		}

		end = std::chrono::steady_clock::now();
		std::cout << "Time for resolvents: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		return resolvents;
	}

	void Model::getEnergies(std::vector<std::vector<double>>& reciever, double_prec direction)
	{
		reciever.reserve(2 * Constants::K_DISCRETIZATION);
		Eigen::SelfAdjointEigenSolver<Matrix_L> solver;
		double_prec k_val = 0;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			k_val = k * L_PI / Constants::K_DISCRETIZATION;
			fillHamiltonian(cos(L_PI * direction) * k_val, sin(L_PI * direction) * k_val);
			solver.compute(hamilton, false);
			reciever.push_back(std::vector<double>(solver.eigenvalues().data(), solver.eigenvalues().data() + solver.eigenvalues().size()));
		}
	}

	void Model::getAllEnergies(std::vector<std::vector<double>>& reciever)
	{
		reciever.resize(4 * Constants::K_DISCRETIZATION, std::vector<double>(2 * Constants::K_DISCRETIZATION));
		Eigen::SelfAdjointEigenSolver<Matrix_L> solver;
		double_prec k_val = 0;
		double_prec l_val = 0;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			k_val = k * L_PI / Constants::K_DISCRETIZATION;
			for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
			{
				l_val = l * L_PI / Constants::K_DISCRETIZATION;
				fillHamiltonian(k_val, l_val);
				solver.compute(hamilton, false);
				reciever[k + Constants::K_DISCRETIZATION][l + Constants::K_DISCRETIZATION] = solver.eigenvalues()(0);

				for (int i = 1; i < 4; i++)
				{
					if (std::abs(solver.eigenvalues()(0) - solver.eigenvalues()(i)) > 1e-8) {
						reciever[k + 3 * Constants::K_DISCRETIZATION][l + Constants::K_DISCRETIZATION] = solver.eigenvalues()(i);
						break;
					}
				}
			}
		}
	}
}