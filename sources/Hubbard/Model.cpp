#include "Model.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <fstream>
#include <thread>
#include <algorithm>
#include <Eigen/SVD>

namespace Hubbard {
	constexpr double SQRT_SALT = 1e-6;
	constexpr double SALT = SQRT_SALT * SQRT_SALT;

	void Model::computeChemicalPotential()
	{
		chemical_potential = U / 2;
	}

	double_prec Model::computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const
	{
		const Eigen::Vector2i l_idx = { x(l), y(l) };
		const Eigen::Vector2i k_idx = { x(k), y(k) };
		Eigen::Vector2i q_idx = { 0,0 };
		std::vector<Eigen::Vector2i> indizes = { l_idx, k_idx, q_idx };
		Eigen::Vector2i momentum_value, coeff_momentum;

		if (term.coefficients.size() > 1) throw std::invalid_argument("Undefined number of coefficients: " + std::to_string(term.coefficients.size()));
		if (term.operators.size() == 0) {
			if (term.coefficients.size() == 1) {
				coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
				return term.multiplicity * computeCoefficient(term.coefficients[0], coeff_momentum);
			}
			return term.multiplicity;
		}

		auto compute_single_sum = [&]() -> double {
			double_prec sumBuffer = 0;
			double_prec returnBuffer = 0;
			for (int q = 0; q < BASIS_SIZE; q++)
			{
				q_idx = { x(q), y(q) };
				indizes = { l_idx, k_idx, q_idx };
				sumBuffer = 1;
				for (size_t i = 0; i < term.operators.size(); i++)
				{
					auto it = wick_map.find(term.operators[i].type);
					if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + term.operators[i].type);
					momentum_value = computeMomentum(term.operators[i].momentum, indizes, { 'l', 'k', 'q' });
					sumBuffer *= expecs[it->second](momentum_value(0), momentum_value(1));
				}
				coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k', 'q' });
				if (term.coefficients.size() == 1) {
					sumBuffer *= computeCoefficient(term.coefficients[0], coeff_momentum);
				}
				returnBuffer += sumBuffer;
			}
			return term.multiplicity * returnBuffer;
		};

		if (term.sum_momenta.size() > 0) {
			if (term.sum_momenta.size() > 1) throw std::invalid_argument("Too many sums: " + term.sum_momenta.size());
			if (term.operators.size() == 1) {
				// bilinear term
				auto it = wick_map.find(term.operators[0].type);
				if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + term.operators[0].type);

				if (term.sum_momenta.size() > 1) throw std::invalid_argument("Term with more than one momentum summation: " + term.sum_momenta.size());
				if (term.delta_momenta.size() == 0) throw std::invalid_argument("There is a summation without delta_kl in a bilinear term.");
				if (term.coefficients.size() == 1) {
					if (term.coefficients.back().momentum.momentum_list.size() == 0) {
						coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
						return term.multiplicity * computeCoefficient(term.coefficients[0], coeff_momentum) * sum_of_all[it->second];
					}
					else {
						return compute_single_sum();
					}
				}
				return term.multiplicity * sum_of_all[it->second];
			}
			if (term.operators.size() == 2) {
				// quartic term
				return compute_single_sum();
			}
			throw std::invalid_argument("There are more than 2 WickOperators: " + term.operators.size());
		}

		double returnBuffer = 1;
		for (size_t i = 0; i < term.operators.size(); i++)
		{
			auto it = wick_map.find(term.operators[i].type);
			if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + term.operators[i].type);
			Eigen::Vector2i momentum_value = computeMomentum(term.operators[i].momentum, indizes, { 'l', 'k' });
			returnBuffer *= expecs[it->second](momentum_value(0), momentum_value(1));
		}
		if (term.coefficients.size() == 1) {
			coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
			return term.multiplicity * computeCoefficient(term.coefficients[0], coeff_momentum) * returnBuffer;
		}
		return term.multiplicity * returnBuffer;
	}

	void Model::fill_M_N()
	{
		M = matrixL::Zero(TOTAL_BASIS, TOTAL_BASIS);
		N = matrixL::Zero(TOTAL_BASIS, TOTAL_BASIS);

#pragma omp parallel for
		for (int i = 0; i < number_of_basis_terms; i++)
		{
			double_prec valueBuffer = 0;
			for (int j = 0; j < number_of_basis_terms; j++)
			{
				// fill N
				for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
					for (int k = 0; k < BASIS_SIZE; k++)
					{
						valueBuffer = 0;
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
							valueBuffer = computeTerm(term, l_buf, k);
							N(j * BASIS_SIZE + l_buf, i * BASIS_SIZE + k) += valueBuffer;
						}
						else {
							for (int l = 0; l < BASIS_SIZE; l++)
							{
								valueBuffer = computeTerm(term, l, k);
								N(j * BASIS_SIZE + l, i * BASIS_SIZE + k) += valueBuffer;
							}
						}
					}
				}

				// fill M
				for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
					for (int k = 0; k < BASIS_SIZE; k++)
					{
						valueBuffer = 0;
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
							valueBuffer = computeTerm(term, l_buf, k);
							M(j * BASIS_SIZE + l_buf, i * BASIS_SIZE + k) += valueBuffer;
						}
						else {
							for (int l = 0; l < BASIS_SIZE; l++)
							{
								valueBuffer = computeTerm(term, l, k);
								M(j * BASIS_SIZE + l, i * BASIS_SIZE + k) += valueBuffer;
							}
						}
					}
				}
			}
		}
	}

	void Model::fill_M_N_xp_basis()
	{
		const std::vector<int> cdw_basis_positions = { 2,3,8,9 };
		const size_t hermitian_offsets[6] = {
			0,									BASIS_SIZE,
			(3 * BASIS_SIZE) / 2,				2 * BASIS_SIZE,
			3 * BASIS_SIZE,						4 * BASIS_SIZE
		};
		const size_t antihermitian_offsets[4] = {
			0,						BASIS_SIZE,
			(3 * BASIS_SIZE) / 2,	2 * BASIS_SIZE
		};

		K_plus = matrixL::Zero(5 * BASIS_SIZE, 5 * BASIS_SIZE);
		K_minus = matrixL::Zero(3 * BASIS_SIZE, 3 * BASIS_SIZE);
		L = matrixL::Zero(5 * BASIS_SIZE, 3 * BASIS_SIZE);

#pragma omp parallel for
		for (int i = 0; i < number_of_basis_terms; i++)
		{
			size_t sum_limit = BASIS_SIZE;
			size_t inner_sum_limit = BASIS_SIZE;
			if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()) {
				inner_sum_limit = BASIS_SIZE;
			}
			else {
				inner_sum_limit = BASIS_SIZE / 2;
			}

			for (int j = 0; j < number_of_basis_terms; j++)
			{
				if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), j) == cdw_basis_positions.end()) {
					sum_limit = BASIS_SIZE;
				}
				else {
					sum_limit = BASIS_SIZE / 2;
				}

				// L
				if (i < 6 && j > 5) {
					for (const auto& term : wicks_N[number_of_basis_terms * i + j]) {
						for (int k = 0; k < sum_limit; k++)
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

								if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()) {
									L(hermitian_offsets[i] + l_buf, antihermitian_offsets[j - 6] + k) += computeTerm(term, l_buf, k);
								}
								else {
									if (l_buf >= BASIS_SIZE / 2) {
										continue;
									}
									L(hermitian_offsets[i] + l_buf, antihermitian_offsets[j - 6] + k) += computeTerm(term, l_buf, k);
								}
							}
							else {
								for (int l = 0; l < inner_sum_limit; l++)
								{
									L(hermitian_offsets[i] + l, antihermitian_offsets[j - 6] + k) += computeTerm(term, l, k);
								}
							}
						} // end k-loop
					} // end term-loop
				}
				
				// K_+ / K_-
				// Ignore the offdiagonal blocks as they are 0
				if (i < 6 && j > 5) continue;
				if (j < 6 && i > 5) continue;
				
				for (const auto& term : wicks_M[number_of_basis_terms * i + j]) {
					for (int k = 0; k < sum_limit; k++)
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

							if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()) {
								if (i < 6) {
									K_plus(hermitian_offsets[i] + l_buf, hermitian_offsets[j] + k) += computeTerm(term, l_buf, k);
								}
								else {
									K_minus(antihermitian_offsets[i - 6] + l_buf, antihermitian_offsets[j - 6] + k) += computeTerm(term, l_buf, k);
								}	
							}
							else {
								if (l_buf >= BASIS_SIZE / 2) {
									continue;
								}
								if (i < 6) {
									K_plus(hermitian_offsets[i] + l_buf, hermitian_offsets[j] + k) += computeTerm(term, l_buf, k);
								}
								else {
									K_minus(antihermitian_offsets[i - 6] + l_buf, antihermitian_offsets[j - 6] + k) += computeTerm(term, l_buf, k);
								}
							}
						}
						else {
							for (int l = 0; l < inner_sum_limit; l++)
							{
								if (i < 6) {
									K_plus(hermitian_offsets[i] + l, hermitian_offsets[j] + k) += computeTerm(term, l, k);
								}
								else {
									K_minus(antihermitian_offsets[i - 6] + l, antihermitian_offsets[j - 6] + k) += computeTerm(term, l, k);
								}
							}
						}
					} // end k-loop
				} // end term-loop
			}
		}
	}

	void Model::initializeParameters()
	{
		this->computeChemicalPotential();
		this->BASIS_SIZE = 4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION;
		if (start_basis_at < 0) {
			// We investigate the special x-p-basis
			this->TOTAL_BASIS = this->BASIS_SIZE * 8;
			this->number_of_basis_terms = 10;
		}
		else {
			this->TOTAL_BASIS = this->BASIS_SIZE * this->number_of_basis_terms;
		}
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
	}

	Model::Model(double _temperature, double _U, int _number_of_basis_terms, int _start_basis_at)
		: temperature(_temperature), U(_U), number_of_basis_terms(_number_of_basis_terms), start_basis_at(_start_basis_at)
	{
		initializeParameters();
	}

	Model::Model(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: temperature(_params.temperature), U(_params.U), number_of_basis_terms(_number_of_basis_terms), start_basis_at(_start_basis_at)
	{
		initializeParameters();
	}

	void Model::loadWick(const std::string& filename)
	{
		wicks_M.resize(number_of_basis_terms * number_of_basis_terms);
		wicks_N.resize(number_of_basis_terms * number_of_basis_terms);
		const int name_offset = (start_basis_at < 0) ? 0 : start_basis_at;

		for (int i = 0; i < number_of_basis_terms; i++)
		{
			for (int j = 0; j < number_of_basis_terms; j++)
			{
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "M_" + std::to_string(j + name_offset) + "_" + std::to_string(i + name_offset) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_M[j * number_of_basis_terms + i].clear();
					ia >> wicks_M[j * number_of_basis_terms + i];
					ifs.close();
				}
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "N_" + std::to_string(j + name_offset) + "_" + std::to_string(i + name_offset) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_N[j * number_of_basis_terms + i].clear();
					ia >> wicks_N[j * number_of_basis_terms + i];
					ifs.close();
				}
			}
		}
	}

	std::unique_ptr<Utility::Resolvent<double_prec>> Model::computeCollectiveModes(std::vector<std::vector<double>>& reciever)
	{
		std::cout << "Gap values:  " << delta_cdw << "  " << delta_sc << "  " << delta_eta << std::endl;
		// First off we need to compute every possible expectation value
		// We use the mean field system's symmetries
		// i.e. there are only the standard SC, CDW, Eta and N operators non-zero
		// spin symmetry and conservation of momentum up to addition of Q holds
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		std::chrono::time_point end = std::chrono::steady_clock::now();

		{
			matrixL rho = matrixL::Zero(4, 4);
			Eigen::SelfAdjointEigenSolver<matrixL> solver;

			expecs = std::vector<matrixL>(4, matrixL::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
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
				}
			}
		}
		computeChemicalPotential();
		std::cout << "Filling of all spin ups = " << sum_of_all[0] / BASIS_SIZE << std::endl;
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
		}
		for (size_t i = 0; i < M.rows(); i++)
		{
			for (size_t j = 0; j < M.cols(); j++)
			{
				if (std::abs(M(i, j)) < 1e-13) {
					M(i, j) = 0;
				}
			}
		}
		for (size_t i = 0; i < N.rows(); i++)
		{
			for (size_t j = 0; j < N.cols(); j++)
			{
				if (std::abs(N(i, j)) < 1e-13) {
					N(i, j) = 0;
				}
			}
		}
		{
#pragma omp parallel for
			for (size_t i = 0; i < M.rows(); i++)
			{
				for (size_t j = i + 1; j < M.cols(); j++)
				{
					if (std::abs(M(i, j) - M(j, i)) > 1e-11) {
						std::cerr << std::scientific << std::setprecision(12) << std::abs(M(i, j) - M(j, i)) << std::endl;
						throw std::invalid_argument("M is not hermitian: " + std::to_string(M(i, j))
							+ " || " + std::to_string(M(j, i)) + "\t\tPosition: " + std::to_string(i) + ", " + std::to_string(j));
					}
					if (std::abs(N(i, j) - N(j, i)) > 1e-11) {
						std::cerr << std::scientific << std::setprecision(12) << std::abs(N(i, j) - N(j, i)) << std::endl;
						throw std::invalid_argument("N is not hermitian: " + std::to_string(N(i, j))
							+ " || " + std::to_string(N(j, i)) + "\t\tPosition: " + std::to_string(i) + ", " + std::to_string(j));
					}
				}
				if (std::abs(N(i, i)) < 1e-8) {
					N(i, i) = 0;
				}
			}
		}

		end = std::chrono::steady_clock::now();
		std::cout << "Time for filling of M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		matrixL buffer = K_plus;
		K_plus = K_minus;
		K_minus = buffer;
		L.transposeInPlace();

		vectorL startingState = vectorL::Zero(K_minus.rows());
		for (size_t i = 0; i < BASIS_SIZE; i++)
		{
			startingState(i) = 1;
		}
		startingState.normalize();
		matrixL N_new;
		// M_new = K_plus
		matrixL solver_matrix;
		matrixL inv_sqrt_mat;

		omp_set_num_threads(8);
		omp_set_nested(1);
		Eigen::initParallel();
		{
#pragma omp parallel sections
			{
#pragma omp section
				{
					std::chrono::time_point begin_in = std::chrono::steady_clock::now();

					Eigen::SelfAdjointEigenSolver<matrixL> solver(K_plus);
					vectorL K_EV = solver.eigenvalues();
					for (size_t i = 0; i < solver.eigenvalues().size(); i++)
					{
						if (K_EV(i) < -SALT) {
							std::cerr << K_EV(i) << std::endl;
							throw std::invalid_argument("K_+ is not positive!  " + std::to_string(K_EV(i)));
						}
						else if (K_EV(i) < SALT) {
							K_EV(i) = SALT;
						}
					}
					K_plus = solver.eigenvectors() * K_EV.asDiagonal() * solver.eigenvectors().adjoint();

					std::chrono::time_point end_in = std::chrono::steady_clock::now();
					std::cout << "Time for getting K_plus: "
						<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
				}
#pragma omp section
				{
					std::chrono::time_point begin_in = std::chrono::steady_clock::now();

					Eigen::SelfAdjointEigenSolver<matrixL> solver(K_minus);
					vectorL k_inv_ev = solver.eigenvalues();
					for (size_t i = 0; i < solver.eigenvalues().size(); i++)
					{
						if (solver.eigenvalues()(i) < -SALT) {
							std::cerr << solver.eigenvalues()(i) << std::endl;
							throw std::invalid_argument("K_- is not positive!  " + std::to_string(solver.eigenvalues()(i)));
						}
						else if (solver.eigenvalues()(i) < SALT) {
							k_inv_ev(i) = (1 / SALT);
						}
						else {
							k_inv_ev(i) = 1. / solver.eigenvalues()(i);
						}
					}
					matrixL buffer_matrix = L * solver.eigenvectors();
					N_new = buffer_matrix * k_inv_ev.asDiagonal() * buffer_matrix.adjoint();

					solver.compute(N_new);
					k_inv_ev = solver.eigenvalues();
					for (size_t i = 0; i < solver.eigenvalues().size(); i++)
					{
						if (solver.eigenvalues()(i) < -SALT) {
							std::cerr << solver.eigenvalues()(i) << std::endl;
							throw std::invalid_argument("N_new is not positive!  " + std::to_string(solver.eigenvalues()(i)));
						}
						else if (k_inv_ev(i) < SALT) {
							k_inv_ev(i) = SQRT_SALT;
						}
						else {
							k_inv_ev(i) = 1. / sqrt(k_inv_ev(i));
						}
					}
					inv_sqrt_mat = solver.eigenvectors() * k_inv_ev.asDiagonal() * solver.eigenvectors().adjoint();
					startingState = inv_sqrt_mat * L * startingState;

					std::chrono::time_point end_in = std::chrono::steady_clock::now();
					std::cout << "Time for getting N_new and its inverse sqrt: "
						<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
				}
			}
		}

		begin = std::chrono::steady_clock::now();
		solver_matrix = inv_sqrt_mat * K_plus * inv_sqrt_mat;
		Eigen::SelfAdjointEigenSolver<matrixL> gen_solver;
		gen_solver.compute(solver_matrix);
		reciever.resize(1);
		{
			vectorL state = gen_solver.eigenvectors().transpose() * startingState;
		
			const double RANGE = 10;
			const int STEPS = 15000;
			for (double z = 0; z < RANGE; z += (RANGE / STEPS))
			{
				std::complex<double_prec> z_tilde = std::complex<double_prec>(z*z, 1e-2);
				Eigen::Vector<std::complex<double_prec>, Eigen::Dynamic> diag = 1. / (z_tilde - gen_solver.eigenvalues().array());
				reciever[0].emplace_back(-(state.dot(diag.asDiagonal() * state)).imag());
			}
		}

		Utility::Resolvent<double_prec> R(startingState);
		R.compute(solver_matrix, 200);

		end = std::chrono::steady_clock::now();
		std::cout << "Time for resolvent: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		return std::make_unique<Utility::Resolvent<double_prec>>(R);
	}

	void Model::getEnergies(std::vector<std::vector<double>>& reciever, double direction)
	{
		reciever.reserve(2 * Constants::K_DISCRETIZATION);
		Eigen::SelfAdjointEigenSolver<matrixL> solver;
		double k_val = 0;
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
		Eigen::SelfAdjointEigenSolver<matrixL> solver;
		double k_val = 0;
		double l_val = 0;
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