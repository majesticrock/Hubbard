#include "SquareGeneral.hpp"

namespace Hubbard::Helper {
	void SquareGeneral::fill_block_M(int i, int j)
	{
		for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0U) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
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
		for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0U) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
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
	void SquareGeneral::setNewModelParameters(Utility::InputFileReader& input, const ModelParameters& modelParameters)
	{
		this->internal_setNewModelParameters(input, modelParameters);
	}
}