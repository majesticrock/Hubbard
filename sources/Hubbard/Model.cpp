#define _USE_MATH_DEFINES
#define _NUM(momentum) (expecs[0](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))
#define _CDW(momentum) (expecs[1](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))
#define _SC(momentum) (expecs[2](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))
#define _ETA(momentum) (expecs[3](x(momentum) + Constants::K_DISCRETIZATION, y(momentum) + Constants::K_DISCRETIZATION))

#include "Model.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <fstream>

namespace Hubbard {
	void Model::computeChemicalPotential()
	{
		chemical_potential = U / 2;
	}

	double Model::computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const
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

		if (term.sum_momenta.size() > 0) {
			if (term.sum_momenta.size() > 1) throw std::invalid_argument("Too many sums: " + term.sum_momenta.size());
			if (term.operators.size() == 1) {
				// bilinear term
				auto it = wick_map.find(term.operators[0].type);
				if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + term.operators[0].type);

				if (term.sum_momenta.size() > 1) throw std::invalid_argument("Term with more than one momentum summation: " + term.sum_momenta.size());
				if (term.delta_momenta.size() == 0) throw std::invalid_argument("There is a summation without delta_kl in a bilinear term.");
				if (term.coefficients.size() == 1) {
					return term.multiplicity * computeCoefficient(term.coefficients[0]) * sum_of_all[it->second];
				}
				return term.multiplicity * sum_of_all[it->second];
			}
			if (term.operators.size() == 2) {
				// quartic term
				double sumBuffer = 0;
				double returnBuffer = 0;
				for (int q = 0; q < BASIS_SIZE; q++)
				{
					q_idx = { x(q), y(q) };
					sumBuffer = 1;
					for (size_t i = 0; i < term.operators.size(); i++)
					{
						auto it = wick_map.find(term.operators[i].type);
						if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + term.operators[i].type);
						indizes = { l_idx, k_idx, q_idx };
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
		M = Eigen::MatrixXd::Zero(TOTAL_BASIS, TOTAL_BASIS);
		N = Eigen::MatrixXd::Zero(TOTAL_BASIS, TOTAL_BASIS);

		double valueBuffer = 0;

		for (int i = 0; i < number_of_basis_terms; i++)
		{
			for (int j = 0; j < number_of_basis_terms; j++)
			{
				// fill N
				for (const auto& term : wicks_N[number_of_basis_terms * i + j]) {
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
							//if (i != j) N(j * BASIS_SIZE + k, i * BASIS_SIZE + l_buf) += valueBuffer;
						}
						else {
							for (int l = 0; l < BASIS_SIZE; l++)
							{
								valueBuffer = computeTerm(term, l, k);
								N(j * BASIS_SIZE + l, i * BASIS_SIZE + k) += valueBuffer;
								//if (i != j) N(j * BASIS_SIZE + k, i * BASIS_SIZE + l) += valueBuffer;
							}
						}
					}
				}

				// fill M
				for (const auto& term : wicks_M[number_of_basis_terms * i + j]) {
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
							M(i * BASIS_SIZE + l_buf, j * BASIS_SIZE + k) += valueBuffer;
							//if (i != j) M(j * BASIS_SIZE + k, i * BASIS_SIZE + l_buf) += valueBuffer;
						}
						else {
							for (int l = 0; l < BASIS_SIZE; l++)
							{
								valueBuffer = computeTerm(term, l, k);
								M(i * BASIS_SIZE + l, j * BASIS_SIZE + k) += valueBuffer;
								//if (i != j) M(j * BASIS_SIZE + k, i * BASIS_SIZE + l) += valueBuffer;
							}
						}
					}
				}
			}
		}

		//return;
		//for (size_t i = 0; i < N.rows(); i++)
		//{
		//	if (std::abs(N(i, i)) < 1e-10) {
		//		size_t numRows = N.rows() - 1;
		//		size_t numCols = N.cols();
		//
		//		N.block(i, 0, numRows - i, numCols) = N.block(i + 1, 0, numRows - i, numCols);
		//		M.block(i, 0, numRows - i, numCols) = M.block(i + 1, 0, numRows - i, numCols);
		//		N.conservativeResize(numRows, numCols);
		//		M.conservativeResize(numRows, numCols);
		//		--numCols;
		//		
		//		N.block(0, i, numRows, numCols - i) = N.block(0, i + 1, numRows, numCols - i);
		//		M.block(0, i, numRows, numCols - i) = M.block(0, i + 1, numRows, numCols - i);
		//		N.conservativeResize(numRows, numCols);
		//		M.conservativeResize(numRows, numCols);
		//	}
		//}
	}

	Model::Model(double _temperature, double _U, int _number_of_basis_terms, int _start_basis_at)
		: temperature(_temperature), U(_U), number_of_basis_terms(_number_of_basis_terms), start_basis_at(_start_basis_at)
	{
		this->chemical_potential = U / 2;
		this->BASIS_SIZE = 4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION;
		this->TOTAL_BASIS = this->BASIS_SIZE * this->number_of_basis_terms;
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
	}

	Model::Model(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: temperature(_params.temperature), U(_params.U), number_of_basis_terms(_number_of_basis_terms), start_basis_at(_start_basis_at)
	{
		this->chemical_potential = U / 2;
		this->BASIS_SIZE = 4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION;
		this->TOTAL_BASIS = this->BASIS_SIZE * this->number_of_basis_terms;
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
	}

	void Model::loadWick(const std::string& filename)
	{
		wicks_M.resize(number_of_basis_terms * number_of_basis_terms);
		wicks_N.resize(number_of_basis_terms * number_of_basis_terms);
		for (int i = 0; i < number_of_basis_terms; i++)
		{
			for (int j = 0; j < number_of_basis_terms; j++)
			{
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "M_" + std::to_string(i + start_basis_at) + "_" + std::to_string(j + start_basis_at) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_M[i * number_of_basis_terms + j].clear();
					ia >> wicks_M[i * number_of_basis_terms + j];
					ifs.close();
				}
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "N_" + std::to_string(i + start_basis_at) + "_" + std::to_string(j + start_basis_at) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_N[i * number_of_basis_terms + j].clear();
					ia >> wicks_N[i * number_of_basis_terms + j];
					ifs.close();
				}
			}
		}
	}

	std::unique_ptr<Utility::Resolvent> Model::computeCollectiveModes(std::vector<std::vector<double>>& reciever)
	{
		std::cout << "Gap values:  " << delta_cdw << "  " << delta_sc << "  " << delta_eta << std::endl;
		//Constants::K_DISCRETIZATION *= 2;
		//BASIS_SIZE *= 4;
		// First off we need to compute every possible expectation value
		// We use the mean field system's symmetries
		// i.e. there are only the standard SC, CDW, Eta and N operators non-zero
		// spin symmetry and conservation of momentum up to addition of Q holds
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		std::chrono::time_point end = std::chrono::steady_clock::now();

		Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

		expecs = std::vector<Eigen::MatrixXd>(4, Eigen::MatrixXd::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
			{
				fillHamiltonian((k * M_PI) / Constants::K_DISCRETIZATION, (l * M_PI) / Constants::K_DISCRETIZATION);
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

		double deviation = 0;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
			{
				double a = unperturbed_energy((k * M_PI) / Constants::K_DISCRETIZATION, (l * M_PI) / Constants::K_DISCRETIZATION) 
					* expecs[1](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION);
				double b = (-U / (2 * BASIS_SIZE)) * sum_of_all[1]
					* (1 - 2 * expecs[0](k + Constants::K_DISCRETIZATION, l + Constants::K_DISCRETIZATION));
				deviation += (a - b) * (a - b);
			}
		}
		std::cout << "Deviation from epsilon-g formula: " << sqrt(deviation) << std::endl;
		
		computeChemicalPotential();
		std::cout << "Filling of all spin ups = " << sum_of_all[0] / BASIS_SIZE << std::endl;
		loadWick("../commutators/wick_");
		end = std::chrono::steady_clock::now();
		std::cout << "Time for expectation values: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		fill_M_N();
		for (size_t i = 0; i < M.rows(); i++)
		{
			for (size_t j = 0; j < M.cols(); j++)
			{
				if (std::abs(M(i, j)) < 1e-13) {
					M(i, j) = 0;
				}
			}
		}

		M += 1e-12 * Eigen::MatrixXd::Identity(M.rows(), M.rows());
		//N += 1e-12 * Eigen::MatrixXd::Identity(M.rows(), M.rows());
		end = std::chrono::steady_clock::now();
		std::cout << "Time for filling of M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		/*begin = std::chrono::steady_clock::now();

		//Eigen::MatrixXd M_T = M.transpose();
		//M += M_T;
		//M *= 0.5;
		const Eigen::MatrixXd test_M = M - M.transpose();
		for (size_t i = 0; i < test_M.rows(); i++)
		{
			for (size_t j = i + 1; j < test_M.cols(); j++)
			{
				if (std::abs(test_M(i, j)) > 1e-6) {
					throw std::invalid_argument("M is not hermitian: " + std::to_string(M(i, j)) + " || " + std::to_string(M(j, i)) + "\t\tPosition: " + std::to_string(i) + ", " + std::to_string(j));
				}
			}
		}
		const Eigen::MatrixXd test_N = N - N.transpose();
		for (size_t i = 0; i < test_N.rows(); i++)
		{
			for (size_t j = i + 1; j < test_N.cols(); j++)
			{
				if (std::abs(test_N(i, j)) > 1e-12) {
					throw std::invalid_argument("N is not hermitian: " + std::to_string(N(i, j)) + " || " + std::to_string(N(j, i)) + "\t\tPosition: " + std::to_string(i) + ", " + std::to_string(j));
				}
			}
			if (std::abs(N(i, i)) < 1e-8) {
				N(i, i) = 0;
			}
		}

		solver.compute(M);
		int singular = 0;
		auto m_ev = solver.eigenvalues();
		for (size_t i = 0; i < solver.eigenvalues().size(); i++)
		{
			if (solver.eigenvalues()(i) < 0) {
				if (solver.eigenvalues()(i) == 0.0 || solver.eigenvalues()(i) == -0.0) continue;
				std::cout << std::scientific << std::setprecision(12) << solver.eigenvalues()(i) << std::endl;
				throw std::invalid_argument(": M is not positive!    " + std::to_string(solver.eigenvalues()(i)));
				
			}
			if (solver.eigenvalues()(i) < 1e-10) {
				if(singular++ == 0) std::cout << "Warning: M is singular! :" << solver.eigenvalues()(i) << std::endl;
				//m_ev(i) = 1e-10;
			}
		}
		M = solver.eigenvectors() * m_ev.asDiagonal() * solver.eigenvectors().transpose();
		std::cout << "Total count of singular eigenvalues of M: " << singular << std::endl;
		end = std::chrono::steady_clock::now();
		std::cout << "Time for checking M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;*/
		begin = std::chrono::steady_clock::now();

		Eigen::LLT<Eigen::MatrixXd> llt_M(M);
		//Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gen_solver;
		//gen_solver.compute(N, M, false);

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> gen_solver;
		Eigen::MatrixXd l_matrix = llt_M.matrixL();
		Eigen::MatrixXd inverse_llt_M = l_matrix.inverse();
		Eigen::MatrixXd solver_matrix = inverse_llt_M * N * inverse_llt_M.transpose();
		gen_solver.compute(solver_matrix, false);

		reciever.resize(1);
		Eigen::VectorXd ev = gen_solver.eigenvalues().real();
		reciever[0] = std::vector<double>(ev.data(), ev.data() + ev.size());

		std::sort(reciever.back().begin(), reciever.back().end());

		end = std::chrono::steady_clock::now();
		std::cout << "Time for solving M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	
		begin = std::chrono::steady_clock::now();
		Eigen::VectorXd startingState = Eigen::VectorXd::Random(M.rows());
		startingState.normalize();
		Utility::Resolvent R(startingState);

		//Eigen::MatrixXd inverse_solve = M.inverse() * N;
		//R.compute(inverse_solve, M, 200);
		R.compute(solver_matrix, Eigen::MatrixXd::Identity(M.rows(), M.cols()), 200);
		R.writeDataToFile("../../data/resolvent.txt");

		end = std::chrono::steady_clock::now();
		std::cout << "Time for resolvent: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		return std::make_unique<Utility::Resolvent>(R);
	}

	void Model::getEnergies(std::vector<std::vector<double>>& reciever, double direction)
	{
		reciever.reserve(2 * Constants::K_DISCRETIZATION);
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
		double k_val = 0;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			k_val = k * M_PI / Constants::K_DISCRETIZATION;
			fillHamiltonian(cos(M_PI * direction) * k_val, sin(M_PI * direction) * k_val);
			solver.compute(hamilton, false);
			reciever.push_back(std::vector<double>(solver.eigenvalues().data(), solver.eigenvalues().data() + solver.eigenvalues().size()));
		}
	}

	void Model::getAllEnergies(std::vector<std::vector<double>>& reciever)
	{
		reciever.resize(4 * Constants::K_DISCRETIZATION, std::vector<double>(2 * Constants::K_DISCRETIZATION));
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
		double k_val = 0;
		double l_val = 0;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			k_val = k * M_PI / Constants::K_DISCRETIZATION;
			for (int l = -Constants::K_DISCRETIZATION; l < Constants::K_DISCRETIZATION; l++)
			{
				l_val = l * M_PI / Constants::K_DISCRETIZATION;
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