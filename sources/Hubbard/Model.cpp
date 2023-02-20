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
#include "../Utility/Resolvent.hpp"

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
		if (term.coefficients.size() > 0) {
			coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
		}

		if (term.coefficients.size() > 1) throw std::invalid_argument("Undefined number of coefficients: " + std::to_string(term.coefficients.size()));
		if (term.operators.size() == 0) {
			if (term.coefficients.size() == 1) {
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
					sumBuffer = 0;
					for (size_t i = 0; i < term.operators.size(); i++)
					{
						auto it = wick_map.find(term.operators[i].type);
						if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + term.operators[i].type);
						indizes = { l_idx, k_idx, q_idx };
						momentum_value = computeMomentum(term.operators[i].momentum, indizes, { 'l', 'k', 'q' });
						sumBuffer += expecs[it->second](momentum_value(0), momentum_value(1));
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

		double returnBuffer = 0;
		for (size_t i = 0; i < term.operators.size(); i++)
		{
			auto it = wick_map.find(term.operators[i].type);
			if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + term.operators[i].type);
			Eigen::Vector2i momentum_value = computeMomentum(term.operators[i].momentum, indizes, { 'l', 'k' });
			returnBuffer += expecs[it->second](momentum_value(0), momentum_value(1));
		}
		if (term.coefficients.size() == 1){
			return term.multiplicity * computeCoefficient(term.coefficients[0], coeff_momentum) * returnBuffer;
		}
		return term.multiplicity * returnBuffer;
	}

	void Model::fill_M_N()
	{
		M = Eigen::MatrixXd::Zero(number_of_basis_terms * BASIS_SIZE, number_of_basis_terms * BASIS_SIZE);
		N = Eigen::MatrixXd::Zero(number_of_basis_terms * BASIS_SIZE, number_of_basis_terms * BASIS_SIZE);

		for (int i = 0; i < number_of_basis_terms; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				// fill N
				for (const auto& term : wicks_N[j * number_of_basis_terms + i]) {
					for (int k = 0; k < BASIS_SIZE; k++)
					{
						if (term.delta_momenta.size() > 0) {
							if (term.delta_momenta.size() > 1) throw std::invalid_argument("Too many deltas: " + term.delta_momenta.size());
							if (term.delta_momenta[0].first.momentum_list.size() != 1) throw std::invalid_argument("First delta list is not of size 1: " + term.delta_momenta[0].first.momentum_list.size());
							if (term.delta_momenta[0].second.momentum_list.size() != 1) throw std::invalid_argument("To be implemented: Second delta list is not of size 1: " + term.delta_momenta[0].second.momentum_list.size());

							int l_buf = k;
							if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
								l_buf += (BASIS_SIZE / 2) + Constants::K_DISCRETIZATION;
							}
							N(k, k) += computeTerm(term, l_buf, k);
						}
						else {
							for (int l = 0; l < BASIS_SIZE; l++)
							{
								N(l, k) += computeTerm(term, l, k);
							}
						}
					}
				}

				// fill M
				for (const auto& term : wicks_M[j * number_of_basis_terms + i]) {
					for (int k = 0; k < BASIS_SIZE; k++)
					{
						if (term.delta_momenta.size() > 0) {
							if (term.delta_momenta.size() > 1) throw std::invalid_argument("Too many deltas: " + term.delta_momenta.size());
							if (term.delta_momenta[0].first.momentum_list.size() != 1) throw std::invalid_argument("First delta list is not of size 1: " + term.delta_momenta[0].first.momentum_list.size());
							if (term.delta_momenta[0].second.momentum_list.size() != 1) throw std::invalid_argument("To be implemented: Second delta list is not of size 1: " + term.delta_momenta[0].second.momentum_list.size());

							int l_buf = k;
							if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
								l_buf += (BASIS_SIZE / 2) + Constants::K_DISCRETIZATION;
							}
							M(k, k) += computeTerm(term, l_buf, k);
						}
						else {
							for (int l = 0; l < BASIS_SIZE; l++)
							{
								M(l, k) += computeTerm(term, l, k);
							}
						}
					}
				}
			}
		}
	}

	Model::Model(double _temperature, double _U)
		: temperature(_temperature), U(_U)
	{
		this->chemical_potential = 0;
		this->BASIS_SIZE = 4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION;
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
	}

	Model::Model(ModelParameters& _params)
		: temperature(_params.temperature), U(_params.U)
	{
		this->chemical_potential = 0;
		this->BASIS_SIZE = 4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION;
		this->delta_cdw = 0.1;
		this->delta_sc = 0.1;
		this->delta_eta = 0.001;
	}

	void Model::loadWick(const std::string& filename)
	{
		wicks_M.resize(number_of_basis_terms * (number_of_basis_terms + 1) / 2);
		wicks_N.resize(number_of_basis_terms * (number_of_basis_terms + 1) / 2);
		for (int i = 0; i < number_of_basis_terms; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "M_" + std::to_string(i) + "_" + std::to_string(j) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_M[j * number_of_basis_terms + i].clear();
					ia >> wicks_M[j * number_of_basis_terms + i];
					ifs.close();
				}
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "N_" + std::to_string(i) + "_" + std::to_string(j) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_N[j * number_of_basis_terms + i].clear();
					ia >> wicks_N[j * number_of_basis_terms + i];
					ifs.close();
				}
			}
		}
	}

	void Model::computeCollectiveModes(std::vector<std::vector<double>>& reciever)
	{
		loadWick("../wick_");
		// First off we need to compute every possible expectation value
		// We use the mean field system's symmetries
		// i.e. there are only the standard SC, CDW, Eta and N operators non-zero
		// spin symmetry and conservation of momentum up to addition of Q holds
		std::chrono::time_point begin = std::chrono::steady_clock::now();
		std::chrono::time_point end = std::chrono::steady_clock::now();

		Eigen::MatrixXd rho = Eigen::MatrixXd::Zero(4, 4);
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

		expecs = std::vector<Eigen::MatrixXd>(4, Eigen::MatrixXd::Zero(2 * Constants::K_DISCRETIZATION, 2 * Constants::K_DISCRETIZATION));
		quartic = std::vector<Eigen::MatrixXd>(4, Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE));

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
		std::cout << sum_of_all[0] / BASIS_SIZE << std::endl;
		computeChemicalPotential();

		end = std::chrono::steady_clock::now();
		std::cout << "Time for expectation values: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		fill_M_N();

		end = std::chrono::steady_clock::now();
		std::cout << "Time for filling of M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		begin = std::chrono::steady_clock::now();

		const Eigen::MatrixXd test = M - M.transpose();
		for (size_t i = 0; i < test.rows(); i++)
		{
			for (size_t j = i + 1; j < test.cols(); j++)
			{
				if (abs(test(i, j)) > 1e-8) {
					std::cout << "Not hermitian: " << test(i, j) << "\t\t" << i << "\t" << j << std::endl;
				}
			}
			if (abs(N(i, i)) < 1e-8) {
				N(i, i) = 0;
			}
		}
		M += 1e-13 * Eigen::MatrixXd::Identity(BASIS_SIZE, BASIS_SIZE);
		//N += 1e-6 * Eigen::MatrixXd::Identity(BASIS_SIZE, BASIS_SIZE);
		Eigen::EigenSolver<Eigen::MatrixXd> gen_solver;
		Eigen::MatrixXd inverse_M = M.completeOrthogonalDecomposition().pseudoInverse();
		Eigen::MatrixXd inverse_N = N.completeOrthogonalDecomposition().pseudoInverse();

		gen_solver.compute(inverse_N * M, false);
		for (size_t i = 0; i < gen_solver.eigenvalues().size(); i++)
		{
			if (abs(gen_solver.eigenvalues()(i).imag()) > 1e-8) {
				std::cerr << "Complex eigenvalue! " << gen_solver.eigenvalues()(i) << std::endl;
			}
		}
		reciever.resize(1);
		Eigen::VectorXd ev = gen_solver.eigenvalues().real();
		reciever[0] = std::vector<double>(ev.data(), ev.data() + ev.size());

		std::sort(reciever.back().begin(), reciever.back().end());

		end = std::chrono::steady_clock::now();
		std::cout << "Time for solving M and N: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

		solver.compute(M);
		for (size_t i = 0; i < solver.eigenvalues().size(); i++)
		{
			if (solver.eigenvalues()(i) < 0) {
				std::cerr << i << ": M is not positive!    " << solver.eigenvalues()(i) << std::endl;
			}
			//if (abs(solver.eigenvalues()(i)) < 1e-5) {
			//	std::cout << solver.eigenvalues()(i) << "\n" << solver.eigenvalues().col(i) << std::endl << std::endl;
			//}
		}
		return;
		begin = std::chrono::steady_clock::now();
		Eigen::VectorXd startingState = Eigen::VectorXd::Ones(BASIS_SIZE);
		startingState.normalize();
		Utility::Resolvent R(startingState);

		Eigen::MatrixXd inverse_solve = M.inverse() * N;
		R.compute(inverse_solve, M, BASIS_SIZE);
		R.writeDataToFile("../../data/resolvent.txt");

		end = std::chrono::steady_clock::now();
		std::cout << "Time for resolvent: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
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
		reciever.resize(2 * Constants::K_DISCRETIZATION, std::vector<double>(2 * Constants::K_DISCRETIZATION));
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
		double k_val = 0;
		double l_val = 0;
		for (int k = -Constants::K_DISCRETIZATION; k < Constants::K_DISCRETIZATION; k++)
		{
			k_val = k * M_PI / Constants::K_DISCRETIZATION;
			for (int l = 0; l < Constants::K_DISCRETIZATION; l++)
			{
				l_val = l * M_PI / Constants::K_DISCRETIZATION;
				fillHamiltonian(k_val, l_val);
				solver.compute(hamilton, false);
				reciever[k + Constants::K_DISCRETIZATION][l] = solver.eigenvalues()(0);

				for (int i = 1; i < 4; i++)
				{
					if (abs(solver.eigenvalues()(0) - solver.eigenvalues()(i)) > 1e-8) {
						reciever[(k + Constants::K_DISCRETIZATION) % (2 * Constants::K_DISCRETIZATION)][l + Constants::K_DISCRETIZATION] = solver.eigenvalues()(i);
						break;
					}
				}
			}
		}
	}
}