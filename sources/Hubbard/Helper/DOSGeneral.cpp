#include "DOSGeneral.hpp"

namespace Hubbard::Helper {
	using DOS = DensityOfStates::Square;

	void DOSGeneral::fillBlock(int i, int j)
	{
		//std::cout << wicks_M[number_of_basis_terms * j + i][0] << std::endl;
		//std::cout << sum_of_all[2].real() << ", " << expecs[2](0, 0) << "\t" << computeTerm(wicks_M[number_of_basis_terms * j + i][0], 0, 0).real() << std::endl;
		//std::cout << sum_of_all

		// fill N
		for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < 2 * DOS::size(); ++k)
			{
				if (term.delta_momenta.size() > 0) {
					int l{};
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						// delta gamma, -gamma'
						l = k + (k < DOS::size() ? DOS::size() : - DOS::size());
					}
					N(j * 2 * DOS::size() + l, i * 2 * DOS::size() + k) += computeTerm(term, l, k) * DOS::values_v(k);
				}
				// melo by to byt 0
			}
		}

		// fill M
		for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
			for (int k = 0; k < 2 * DOS::size(); ++k)
			{
				if (term.delta_momenta.size() > 0) {
					int l{};
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						// delta gamma, -gamma'
						l = k + (k < DOS::size() ? DOS::size() : -DOS::size());
					}
					M(j * 2 * DOS::size() + l, i * 2 * DOS::size() + k) += computeTerm(term, l, k) * DOS::values_v(k);
				}
				// melo by to byt 0
			}
		}
	}
}