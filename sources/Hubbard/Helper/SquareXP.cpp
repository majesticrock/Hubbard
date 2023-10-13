#include "SquareXP.hpp"
#include <array>

namespace Hubbard::Helper {
	void SquareXP::fillBlock(int i, int j)
	{
		constexpr std::array<int, 4> cdw_basis_positions{ 2,3,8,9 };
		const int hermitian_offsets[6] = {
			0,									Constants::BASIS_SIZE,
			(3 * Constants::BASIS_SIZE) / 2,	2 * Constants::BASIS_SIZE,
			3 * Constants::BASIS_SIZE,			4 * Constants::BASIS_SIZE
		};
		const int antihermitian_offsets[4] = {
			0,									Constants::BASIS_SIZE,
			(3 * Constants::BASIS_SIZE) / 2,	2 * Constants::BASIS_SIZE
		};

		const int sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), j) == cdw_basis_positions.end()
			? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;
		const int inner_sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()
			? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;

		// L
		if (i < 6 && j > 5) {
			for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
				for (int k = 0; k < sum_limit; ++k)
				{
					if (term.delta_momenta.size() > 0) {
						int l{ k };
						if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
							l = addQTo(k);
						}

						if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()) {
							L(hermitian_offsets[i] + k, antihermitian_offsets[j - 6] + l) += computeRealTerm(term, k, l);
						}
						else {
							if (l >= Constants::BASIS_SIZE / 2) {
								continue;
							}
							L(hermitian_offsets[i] + k, antihermitian_offsets[j - 6] + l) += computeRealTerm(term, k, l);
						}
					}
					else {
						for (int l = 0; l < inner_sum_limit; l++)
						{
							L(hermitian_offsets[i] + k, antihermitian_offsets[j - 6] + l) += computeRealTerm(term, k, l);
						}
					}
				} // end k-loop
			} // end term-loop
		}

		// K_+ / K_-
		// Ignore the offdiagonal blocks as they are 0
		if (i < 6 && j > 5) return;
		if (j < 6 && i > 5) return;

		for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
			for (int k = 0; k < sum_limit; ++k)
			{
				if (term.delta_momenta.size() > 0U) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
					}

					if (std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()) {
						if (i < 6) {
							K_plus(hermitian_offsets[i] + k, hermitian_offsets[j] + l) += computeRealTerm(term, k, l);
						}
						else {
							K_minus(antihermitian_offsets[i - 6] + k, antihermitian_offsets[j - 6] + l) += computeRealTerm(term, k, l);
						}
					}
					else {
						if (l >= Constants::BASIS_SIZE / 2) {
							continue;
						}
						if (i < 6) {
							K_plus(hermitian_offsets[i] + k, hermitian_offsets[j] + l) += computeRealTerm(term, k, l);
						}
						else {
							K_minus(antihermitian_offsets[i - 6] + k, antihermitian_offsets[j - 6] + l) += computeRealTerm(term, k, l);
						}
					}
				}
				else {
					for (int l = 0; l < inner_sum_limit; ++l)
					{
						if (i < 6) {
							K_plus(hermitian_offsets[i] + k, hermitian_offsets[j] + l) += computeRealTerm(term, k, l);
						}
						else {
							K_minus(antihermitian_offsets[i - 6] + k, antihermitian_offsets[j - 6] + l) += computeRealTerm(term, k, l);
						}
					}
				}
			} // end k-loop
		} // end term-loop
	}
}