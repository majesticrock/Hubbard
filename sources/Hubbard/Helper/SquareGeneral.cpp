#include "SquareGeneral.hpp"

namespace Hubbard::Helper {
	void SquareGeneral::fillBlock(int i, int j)
	{
		// fill N
		for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0) {
					int l_buf = k;
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						Eigen::Vector2i l_buf_vec = { x(k), y(k) };
						l_buf_vec(0) += Constants::K_DISCRETIZATION;
						l_buf_vec(1) += Constants::K_DISCRETIZATION;
						clean_factor_2pi(l_buf_vec);
						l_buf = l_buf_vec(0) * 2 * Constants::K_DISCRETIZATION + l_buf_vec(1);
					}
					N(j * Constants::BASIS_SIZE + l_buf, i * Constants::BASIS_SIZE + k) += computeTerm(term, l_buf, k);
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; l++)
					{
						N(j * Constants::BASIS_SIZE + l, i * Constants::BASIS_SIZE + k) += computeTerm(term, l, k);
					}
				}
			}
		}
		
		// fill M
		for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
			//if (term.delta_momenta.size() > 0 && term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
			//	static int counter;
			//	int k = 0;
			//	int l_buf = addQTo(k);
			//	if (term.delta_momenta.size() > 0) {
			//		auto buf = computeTerm(term, k, l_buf).real();
			//		if(std::abs(buf) > 1e-12){
			//			std::cout << "$" << counter++ << " || " << term << " = " << buf << "$\n" << std::endl;
			//		}
			//	}
			//}
			
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0) {
					int l_buf = k;
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l_buf = addQTo(k);
					}
					M(j * Constants::BASIS_SIZE + l_buf, i * Constants::BASIS_SIZE + k) += computeTerm(term, l_buf, k);
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; l++)
					{
						M(j * Constants::BASIS_SIZE + l, i * Constants::BASIS_SIZE + k) += computeTerm(term, l, k);
					}
				}
			}
		}
	}
}