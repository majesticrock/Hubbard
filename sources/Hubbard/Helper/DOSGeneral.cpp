#include "DOSGeneral.hpp"

namespace Hubbard::Helper {
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
						l = this->model->shiftByQ(k);
					}
					N(j + k * number_of_basis_terms, i + l * number_of_basis_terms) += computeTerm(term, l, k) * this->approximate_dos[k]* this->approximate_dos[k];
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; ++l)
					{
						N(j + k * number_of_basis_terms, i + l * number_of_basis_terms) += computeTerm(term, l, k) * this->approximate_dos[k] * this->approximate_dos[l];
					}
				}
			}
		}

		// fill M
		for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
			//if (term.delta_momenta.size() > 0 && term.delta_momenta[0].first.add_Q == term.delta_momenta[0].second.add_Q) {
			//	int pos = 0;
			//	auto buf = computeTerm(term, pos, pos + Constants::HALF_BASIS).real();
			//	static int counter;
			//	if(std::abs(buf) > 1e-12){
			//		std::cout << "$" << counter++ << " " << DOS::abscissa_v(pos) << " || " << term << " = " << buf << "$\n" << std::endl;
			//	}
			//}
			for (int k = 0; k < Constants::BASIS_SIZE; ++k)
			{
				if (term.delta_momenta.size() > 0) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						// delta gamma, -gamma'
						l = this->model->shiftByQ(k);
					}
					M(j + k * number_of_basis_terms, i + l * number_of_basis_terms) += computeTerm(term, l, k) * this->approximate_dos[k] * this->approximate_dos[k];
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; ++l)
					{
						M(j + k * number_of_basis_terms, i + l * number_of_basis_terms) += computeTerm(term, l, k) * this->approximate_dos[k] * this->approximate_dos[l];
					}
				}
			}
		}
	}
}