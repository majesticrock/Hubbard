#include "SquareGeneral.hpp"

namespace Hubbard::Helper {
	void SquareGeneral::fillBlock(int i, int j)
	{
		// fill N
		for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0) {
					if (term.delta_momenta.size() > 1) throw std::invalid_argument("Too many deltas: " + term.delta_momenta.size());
					if (term.delta_momenta[0].first.momentum_list.size() != 1) throw std::invalid_argument("First delta list is not of size 1: " + term.delta_momenta[0].first.momentum_list.size());
					if (term.delta_momenta[0].second.momentum_list.size() != 1) throw std::invalid_argument("To be implemented: Second delta list is not of size 1: " + term.delta_momenta[0].second.momentum_list.size());

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
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0) {
					if (term.delta_momenta.size() > 1) throw std::invalid_argument("Too many deltas: " + term.delta_momenta.size());
					if (term.delta_momenta[0].first.momentum_list.size() != 1) throw std::invalid_argument("First delta list is not of size 1: " + term.delta_momenta[0].first.momentum_list.size());
					if (term.delta_momenta[0].second.momentum_list.size() != 1) throw std::invalid_argument("To be implemented: Second delta list is not of size 1: " + term.delta_momenta[0].second.momentum_list.size());

					int l_buf = k;
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						Eigen::Vector2i l_buf_vec = { x(k), y(k) };
						l_buf_vec(0) += Constants::K_DISCRETIZATION;
						l_buf_vec(1) += Constants::K_DISCRETIZATION;
						clean_factor_2pi(l_buf_vec);
						l_buf = l_buf_vec(0) * 2 * Constants::K_DISCRETIZATION + l_buf_vec(1);
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