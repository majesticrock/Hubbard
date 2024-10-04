#include "SquareXP.hpp"

namespace Hubbard::Helper {
	// K_+ / K_-
	void SquareXP::fill_block_M(int i, int j)
	{
		const int sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()
			? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;
		const int inner_sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), j) == cdw_basis_positions.end()
			? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;

		for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
			for (int k = 0; k < sum_limit; ++k)
			{
				if (term.delta_momenta.size() > 0U) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
					}
					if (term.delta_momenta[0].second.momentum_list.size() == 2U) {
						l = add_index_momenta(l, this->mode_momentum, -term.delta_momenta[0].second.momentum_list[1].first);
					}
					if (l >= inner_sum_limit) {
						continue;
					}

					if (i < hermitian_size) {
						K_plus(hermitian_offsets[i] + k, hermitian_offsets[j] + l) += computeRealTerm(term, k, l);
					}
					else {
						K_minus(antihermitian_offsets[i - hermitian_size] + k, antihermitian_offsets[j - hermitian_size] + l) += computeRealTerm(term, k, l);
					}
				}
				else {
					for (int l = 0; l < inner_sum_limit; ++l)
					{
						if (i < hermitian_size) {
							K_plus(hermitian_offsets[i] + k, hermitian_offsets[j] + l) += computeRealTerm(term, k, l);
						}
						else {
							K_minus(antihermitian_offsets[i - hermitian_size] + k, antihermitian_offsets[j - hermitian_size] + l) += computeRealTerm(term, k, l);
						}
					}
				}
			} // end k-loop
		} // end term-loop
	}

	// L
	void SquareXP::fill_block_N(int i, int j)
	{
		const int sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()
			? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;
		const int inner_sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), j) == cdw_basis_positions.end()
			? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;

		for (const auto& term : wicks.N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < sum_limit; ++k)
			{
				if (term.delta_momenta.size() > 0U) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
					}
					if (term.delta_momenta[0].second.momentum_list.size() == 2U) {
						l = add_index_momenta(l, this->mode_momentum, -term.delta_momenta[0].second.momentum_list[1].first);
					}
					if (l >= inner_sum_limit) {
						continue;
					}

					L(hermitian_offsets[i] + k, antihermitian_offsets[j - hermitian_size] + l) += computeRealTerm(term, k, l);
				}
				else {
					for (int l = 0; l < inner_sum_limit; l++)
					{
						L(hermitian_offsets[i] + k, antihermitian_offsets[j - hermitian_size] + l) += computeRealTerm(term, k, l);
					}
				}
			} // end k-loop
		} // end term-loop
	}
	void SquareXP::setNewModelParameters(Utility::InputFileReader& input, const Models::ModelParameters& modelParameters)
	{
		this->internal_setNewModelParameters(input, modelParameters);
	}
}