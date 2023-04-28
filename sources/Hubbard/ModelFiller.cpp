#include "Model.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <fstream>
#include <algorithm>

namespace Hubbard{
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

		auto compute_single_sum = [&]() -> double_prec {
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

		double_prec returnBuffer = 1;
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
		M = Matrix_L::Zero(TOTAL_BASIS, TOTAL_BASIS);
		N = Matrix_L::Zero(TOTAL_BASIS, TOTAL_BASIS);

#pragma omp parallel for
		for (int i = 0; i < number_of_basis_terms; i++)
		{
			for (int j = 0; j < number_of_basis_terms; j++)
			{
				// fill N
				for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
					for (int k = 0; k < BASIS_SIZE; k++)
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
							N(j * BASIS_SIZE + l_buf, i * BASIS_SIZE + k) += computeTerm(term, l_buf, k);
						}
						else {
							for (int l = 0; l < BASIS_SIZE; l++)
							{
								N(j * BASIS_SIZE + l, i * BASIS_SIZE + k) += computeTerm(term, l, k);
							}
						}
					}
				}

				// fill M
				for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
					for (int k = 0; k < BASIS_SIZE; k++)
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
							M(j * BASIS_SIZE + l_buf, i * BASIS_SIZE + k) += computeTerm(term, l_buf, k);
						}
						else {
							for (int l = 0; l < BASIS_SIZE; l++)
							{
								M(j * BASIS_SIZE + l, i * BASIS_SIZE + k) += computeTerm(term, l, k);
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

		K_plus = Matrix_L::Zero(5 * BASIS_SIZE, 5 * BASIS_SIZE);
		K_minus = Matrix_L::Zero(3 * BASIS_SIZE, 3 * BASIS_SIZE);
		L = Matrix_L::Zero(5 * BASIS_SIZE, 3 * BASIS_SIZE);

		//#pragma omp parallel for
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
}