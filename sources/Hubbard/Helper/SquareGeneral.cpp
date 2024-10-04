#include "SquareGeneral.hpp"

namespace Hubbard::Helper {
	void SquareGeneral::fill_block_M(int i, int j)
	{
		for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0U) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
					}
					if (term.delta_momenta[0].second.momentum_list.size() == 2U) {
						l = add_index_momenta(l, this->mode_momentum, -term.delta_momenta[0].second.momentum_list[1].first);
					}
					M(i + k * number_of_basis_terms, j + l * number_of_basis_terms) += computeTerm(term, k, l);
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; l++)
					{
						M(i + k * number_of_basis_terms, j + l * number_of_basis_terms) += computeTerm(term, k, l);
					}
				}
			}
		}
	}

	void SquareGeneral::fill_block_N(int i, int j)
	{
		for (const auto& term : wicks.N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0U) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
					}
					if (term.delta_momenta[0].second.momentum_list.size() == 2U) {
						l = add_index_momenta(l, this->mode_momentum, -term.delta_momenta[0].second.momentum_list[1].first);
					}
					N(i + k * number_of_basis_terms, j + l * number_of_basis_terms) += computeTerm(term, k, l);
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; l++)
					{
						N(i + k * number_of_basis_terms, j + l * number_of_basis_terms) += computeTerm(term, k, l);
					}
				}
			}
		}
	}
	void SquareGeneral::setNewModelParameters(Utility::InputFileReader& input, const Models::ModelParameters& modelParameters)
	{
		this->internal_setNewModelParameters(input, modelParameters);
	}
}