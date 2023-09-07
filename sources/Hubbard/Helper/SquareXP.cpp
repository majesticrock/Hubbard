#include "SquareXP.hpp"

namespace Hubbard::Helper {
	void SquareXP::fillBlock(int i, int j)
	{
		const std::vector<int> cdw_basis_positions = { 2,3,8,9 };
		const int hermitian_offsets[6] = {
			0,									Constants::BASIS_SIZE,
			(3 * Constants::BASIS_SIZE) / 2,				2 * Constants::BASIS_SIZE,
			3 * Constants::BASIS_SIZE,						4 * Constants::BASIS_SIZE
		};
		const int antihermitian_offsets[4] = {
			0,						Constants::BASIS_SIZE,
			(3 * Constants::BASIS_SIZE) / 2,	2 * Constants::BASIS_SIZE
		};

		int sum_limit = Constants::BASIS_SIZE;
		int inner_sum_limit = Constants::BASIS_SIZE;
		if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()) {
			inner_sum_limit = Constants::BASIS_SIZE;
		}
		else {
			inner_sum_limit = Constants::BASIS_SIZE / 2;
		}
		if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), j) == cdw_basis_positions.end()) {
			sum_limit = Constants::BASIS_SIZE;
		}
		else {
			sum_limit = Constants::BASIS_SIZE / 2;
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
							L(hermitian_offsets[i] + l_buf, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l_buf, k);
						}
						else {
							if (l_buf >= Constants::BASIS_SIZE / 2) {
								continue;
							}
							L(hermitian_offsets[i] + l_buf, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l_buf, k);
						}
					}
					else {
						for (int l = 0; l < inner_sum_limit; l++)
						{
							L(hermitian_offsets[i] + l, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l, k);
						}
					}
				} // end k-loop
			} // end term-loop
		}

		// K_+ / K_-
		// Ignore the offdiagonal blocks as they are 0
		if (i < 6 && j > 5) return;
		if (j < 6 && i > 5) return;

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
							K_plus(hermitian_offsets[i] + l_buf, hermitian_offsets[j] + k) += computeRealTerm(term, l_buf, k);
						}
						else {
							K_minus(antihermitian_offsets[i - 6] + l_buf, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l_buf, k);
						}
					}
					else {
						if (l_buf >= Constants::BASIS_SIZE / 2) {
							continue;
						}
						if (i < 6) {
							K_plus(hermitian_offsets[i] + l_buf, hermitian_offsets[j] + k) += computeRealTerm(term, l_buf, k);
						}
						else {
							K_minus(antihermitian_offsets[i - 6] + l_buf, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l_buf, k);
						}
					}
				}
				else {
					for (int l = 0; l < inner_sum_limit; l++)
					{
						if (i < 6) {
							K_plus(hermitian_offsets[i] + l, hermitian_offsets[j] + k) += computeRealTerm(term, l, k);
						}
						else {
							K_minus(antihermitian_offsets[i - 6] + l, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l, k);
						}
					}
				}
			} // end k-loop
		} // end term-loop
	}
}