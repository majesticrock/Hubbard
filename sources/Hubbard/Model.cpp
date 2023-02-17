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

	void Model::compute_quartics()
	{
		const Eigen::Vector2i vec_Q = { Constants::K_DISCRETIZATION, Constants::K_DISCRETIZATION };

#pragma omp parallel for
		for (int k = 0; k < BASIS_SIZE; k++)
		{
			Eigen::Vector2i vec_k, vec_l;
			vec_k = { x(k), y(k) };
			for (int l = 0; l < BASIS_SIZE; l++)
			{
				vec_l = { x(l), y(l) };
				quartic[0](k, l) += sc_type(vec_l, -vec_l) * sc_type(-vec_k, vec_k);
				quartic[0](k, l) += sc_type(vec_l + vec_Q, -vec_l) * sc_type(-vec_k, vec_k + vec_Q);
				quartic[0](k, l) -= cdw_type(-vec_k, -vec_k) * cdw_type(-vec_l, -vec_l);
				quartic[0](k, l) -= cdw_type(-vec_k + vec_Q, -vec_k) * cdw_type(-vec_l, -vec_l + vec_Q);
				if (k == l) { // Factor 2 because of the spin sum
					quartic[0](k, l) += 2 * sum_of_all[0] * cdw_type(-vec_l, -vec_k);
				}
				else if (vec_k == vec_l + vec_Q || vec_k == vec_l - vec_Q) {
					quartic[0](k, l) += 2 * sum_of_all[1] * cdw_type(-vec_l, -vec_k);;
				}

				quartic[1](k, l) += sc_type(vec_l, -vec_l) * sc_type(-vec_k, vec_k);
				quartic[1](k, l) += sc_type(vec_l, -vec_l - vec_Q) * sc_type(-vec_k, vec_k - vec_Q);
				quartic[1](k, l) += cdw_type(vec_l, vec_l) * cdw_type(-vec_k, -vec_k);
				quartic[1](k, l) += cdw_type(vec_l, vec_l - vec_Q) * cdw_type(-vec_k - vec_Q, -vec_k);

				quartic[3](k, l) += cdw_type(-vec_k, -vec_k) * cdw_type(-vec_l, -vec_l);
				quartic[3](k, l) += cdw_type(-vec_k - vec_Q, -vec_k) * cdw_type(-vec_l, -vec_l - vec_Q);
				if (k == l) {
					quartic[3](k, l) -= sum_of_all[0] * cdw_type(-vec_l, -vec_k);
				}
				else if (vec_k == vec_l + vec_Q || vec_k == vec_l - vec_Q) {
					quartic[3](k, l) -= sum_of_all[1] * cdw_type(-vec_l, -vec_k);
				}
			}
			// the formulas include delta_kl, hence we dont need to recompute the same term for each l
			quartic[2](k, k) += sc_type(vec_k, -vec_k) * sum_of_all[2];
			quartic[2](k, k) += sc_type(vec_k, -vec_k + vec_Q) * sum_of_all[3];
			quartic[2](k, k) += cdw_type(vec_k, vec_k) * sum_of_all[0];
			quartic[2](k, k) += cdw_type(vec_k + vec_Q, vec_k) * sum_of_all[1];
		}
	}

	double Model::computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const
	{
		if (term.operators.size() == 0) return term.multiplicity;
		const Eigen::Vector2i l_idx = { x(l) + Constants::K_DISCRETIZATION, y(l) + Constants::K_DISCRETIZATION };
		const Eigen::Vector2i k_idx = { x(k) + Constants::K_DISCRETIZATION, y(k) + Constants::K_DISCRETIZATION };

		auto computeMomentum_kl = [&](const SymbolicOperators::WickOperator& op) -> Eigen::Vector2i {
			Eigen::Vector2i buffer = { 0,0 };
			int idx = op.momentum.isUsed('k');
			int jdx = op.momentum.isUsed('l');

			if (idx > -1) {
				buffer += op.momentum.momentum_list[idx].first * k_idx;
			}
			if (jdx > -1) {
				buffer += op.momentum.momentum_list[jdx].first * l_idx;
			}
			if (op.momentum.add_Q) {
				buffer(0) += Constants::K_DISCRETIZATION;
				buffer(1) += Constants::K_DISCRETIZATION;
			}
			clean_factor_2pi(buffer);
			return buffer;
		};

		if (term.sum_momenta.size() > 0) {
			if (term.operators.size() == 1) {
				// bilinear term
				auto it = wick_map.find(term.operators[0].type);
				if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + term.operators[0].type);

				if (term.sum_momenta.size() > 1) throw std::invalid_argument("Term with more than one momentum summation: " + term.sum_momenta.size());
				if (term.delta_momenta.size() == 0) throw std::invalid_argument("There is a summation without delta_kl in a bilinear term.");
				return term.multiplicity * sum_of_all[it->second];
			}
			if (term.operators.size() == 2) {
				// quartic term

			}
			throw std::invalid_argument("There are more than 2 WickOperators: " + term.operators.size());
		}

		double returnBuffer = 0;
		for (size_t i = 0; i < term.operators.size(); i++)
		{
			auto it = wick_map.find(term.operators[i].type);
			if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + term.operators[i].type);
			Eigen::Vector2i momentum_value = computeMomentum_kl(term.operators[i]);
			returnBuffer += expecs[it->second](momentum_value(0), momentum_value(1));
		}
		return term.multiplicity * returnBuffer;
	}

	void Model::fill_M_N()
	{
		M = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);
		N = Eigen::MatrixXd::Zero(BASIS_SIZE, BASIS_SIZE);

		for (int k = 0; k < BASIS_SIZE; k++)
		{
			N(k, k) = -(1 - 2 * _NUM(k));

			M(k, k) += 2 * (unperturbed_energy((M_PI * x(k)) / Constants::K_DISCRETIZATION, (M_PI * y(k)) / Constants::K_DISCRETIZATION)
				- chemical_potential) * (1 - 2 * _NUM(k));

			M(k, k) += (2 * U / BASIS_SIZE) * (sum_of_all[0] - quartic[2](k, k));

			for (int l = 0; l < BASIS_SIZE; l++)
			{
				M(l, k) += (U / BASIS_SIZE) * (1 - 2 * (_NUM(l) + _NUM(k) - quartic[1](l, k) - quartic[3](l, k)));
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

	void Model::computeCollectiveModes(std::vector<std::vector<double>>& reciever)
	{
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

		compute_quartics();

		end = std::chrono::steady_clock::now();
		std::cout << "Time for quartic values: "
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
		R.writeDataToFile("../../dataresolvent.txt");

		end = std::chrono::steady_clock::now();
		std::cout << "Time for resolvent: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	}

	void Model::loadWick(const std::string& filename)
	{
		// create an input file stream and a text archive to deserialize the vector
		std::ifstream ifs(filename);
		boost::archive::text_iarchive ia(ifs);
		wicks.clear();
		ia >> wicks;
		ifs.close();
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

	Model::ModelParameters::ModelParameters(double _temperature, double _U, double _V, double _global_step, double _second_step,
		std::string _global_iterator_type, std::string _second_iterator_type)
		: global_iterator_type(_global_iterator_type), second_iterator_type(_second_iterator_type),
		global_step(_global_step), second_step(_second_step), temperature(_temperature), U(_U), V(_V)
	{
		if (second_iterator_type == "T") {
			second_it_min = temperature;
		}
		else if (second_iterator_type == "U") {
			second_it_min = U;
		}
		else if (second_iterator_type == "V") {
			second_it_min = V;
		}
		else {
			second_it_min = 0;
		}
	}

	// // // // // // // // // // // //
	void Model::ModelParameters::incrementer(std::string& s, const double step)
	{
		if (s == "T") {
			temperature += step;
		}
		else if (s == "U") {
			U += step;
		}
		else if (s == "V") {
			V += step;
		}
	}
	void Model::ModelParameters::setSecondIterator(int it_num)
	{
		if (second_iterator_type == "T") {
			temperature = second_it_min + it_num * second_step;
		}
		else if (second_iterator_type == "U") {
			U = second_it_min + it_num * second_step;
		}
		else if (second_iterator_type == "V") {
			V = second_it_min + it_num * second_step;
		}
	}
	void Model::ModelParameters::incrementGlobalIterator()
	{
		incrementer(global_iterator_type, global_step);
		setSecondIterator(0);
	}
	void Model::ModelParameters::incrementSecondIterator()
	{
		incrementer(second_iterator_type, second_step);
	}
	void Model::ModelParameters::printGlobal() const
	{
		std::cout << global_iterator_type << " = " << getGlobal();
	}
	void Model::data_set::print() const {
		std::cout << delta_cdw << "\t" << delta_sc << "\t" << delta_eta
			<< "\t" << sqrt(delta_cdw * delta_cdw + delta_sc * delta_sc + delta_eta * delta_eta) << std::endl;
	}
}