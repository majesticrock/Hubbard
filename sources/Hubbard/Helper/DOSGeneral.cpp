#include "DOSGeneral.hpp"

namespace Hubbard::Helper {
	using DOS = DensityOfStates::Square;

	void DOSGeneral::fillBlock(int i, int j)
	{
		// fill N
		for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < DOS::size(); ++k)
			{
				if (term.delta_momenta.size() > 0) {
					int l{};
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						// delta gamma, -gamma'
						l = k + DOS::size();
					}
					N(j * number_of_basis_terms + l, i * number_of_basis_terms + k) += computeTerm(term, l, k);
				}
				// melo by to byt 0
			}
		}

		// fill M
		for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
			for (int k = 0; k < DOS::size(); ++k)
			{
				if (term.delta_momenta.size() > 0) {
					int l{};
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						// delta gamma, -gamma'
						l = k + DOS::size();
					}
					N(j * number_of_basis_terms + l, i * number_of_basis_terms + k) += computeTerm(term, l, k);
				}
				// melo by to byt 0
			}
		}
	}
}