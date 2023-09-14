#include "SquareGeneral.hpp"

namespace Hubbard::Helper {
	void SquareGeneral::fillBlock(int i, int j)
	{
		// fill N
		for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
					}
					N(j + l * number_of_basis_terms, i + k * number_of_basis_terms) += computeTerm(term, l, k);
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; l++)
					{
						N(j + l * number_of_basis_terms, i + k * number_of_basis_terms) += computeTerm(term, l, k);
					}
				}
			}
		}
		
		// fill M
		for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
			//if (term.delta_momenta.size() > 0 && term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
			//	static int counter;
			//	int k = 31;
			//	int l = addQTo(k);
			//	if (term.delta_momenta.size() > 0) {
			//		auto buf = computeTerm(term, k, l).real();
			//		if(std::abs(buf) > 1e-12){
			//			std::cout << "$" << counter++ << " " << gammaFromIndex(k) << " || " << term << " = " << buf << "$\n" << std::endl;
			//		}
			//	}
			//}
			
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
					}
					M(j + l * number_of_basis_terms, i + k * number_of_basis_terms) += computeTerm(term, l, k);
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; l++)
					{
						M(j + l * number_of_basis_terms, i + k * number_of_basis_terms) += computeTerm(term, l, k);
					}
				}
			}
		}
	}
}