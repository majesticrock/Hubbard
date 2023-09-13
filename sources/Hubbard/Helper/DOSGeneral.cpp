#include "DOSGeneral.hpp"

namespace Hubbard::Helper {
	using DOS = DensityOfStates::Square;
	void DOSGeneral::fillBlock(int i, int j)
	{
		// fill N
		for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; ++k)
			{
				if (term.delta_momenta.size() > 0) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						// delta gamma, -gamma'
						l += (k < DOS::size() ? DOS::size() : -DOS::size());
					}
					N(j * Constants::BASIS_SIZE + l, i * Constants::BASIS_SIZE + k) += computeTerm(term, l, k) * DOS::values_v(k);
				}
				// melo by to byt 0
			}
		}

		// fill M
		for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
			//if (term.delta_momenta.size() > 0 && term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
			//	auto buf = computeTerm(term, DOS::size() - 1, 2 * DOS::size() - 1).real();
			//	static int counter;
			//	if(std::abs(buf) > 1e-12){
			//		std::cout << "$" << counter++ << " || " << term << " = " << buf << "$\n" << std::endl;
			//	}
			//}

			for (int k = 0; k < Constants::BASIS_SIZE; ++k)
			{
				if (term.delta_momenta.size() > 0) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						// delta gamma, -gamma'
						l += (k < DOS::size() ? DOS::size() : -DOS::size());
					}
					M(j * Constants::BASIS_SIZE + l, i * Constants::BASIS_SIZE + k) += computeTerm(term, l, k) * DOS::values_v(k);
				}
				// melo by to byt 0
			}
		}
	}
}