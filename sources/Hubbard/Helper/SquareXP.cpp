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
						int l{ k };
						if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
							l = addQTo(k);
						}

						if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()) {
							L(hermitian_offsets[i] + l, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l, k);
						}
						else {
							if (l >= Constants::BASIS_SIZE / 2) {
								continue;
							}
							L(hermitian_offsets[i] + l, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l, k);
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
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
					}

					if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()) {
						if (i < 6) {
							K_plus(hermitian_offsets[i] + l, hermitian_offsets[j] + k) += computeRealTerm(term, l, k);
						}
						else {
							K_minus(antihermitian_offsets[i - 6] + l, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l, k);
						}
					}
					else {
						if (l >= Constants::BASIS_SIZE / 2) {
							continue;
						}
						if (i < 6) {
							K_plus(hermitian_offsets[i] + l, hermitian_offsets[j] + k) += computeRealTerm(term, l, k);
						}
						else {
							K_minus(antihermitian_offsets[i - 6] + l, antihermitian_offsets[j - 6] + k) += computeRealTerm(term, l, k);
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