#include "SquareGeneral.hpp"

namespace Hubbard::Helper {
	enum BrokenSymmetry : int { None = 0, Momentum = 1, ParticleNumber = 2, Spin = 4 };
	constexpr BrokenSymmetry breaks_symmetry(int i, int j) {
		int sym{};
		if(i == 0 || i == 1 || i == 6 || i == 7) {
			sym ^= BrokenSymmetry::ParticleNumber;
		}
		if(i >= 4 && i <= 7) {
			sym ^= BrokenSymmetry::Momentum;
		}
		if(j == 0 || j == 1 || j == 6 || j == 7) {
			sym ^= BrokenSymmetry::ParticleNumber;
		}
		if(j >= 4 && j <= 7) {
			sym ^= BrokenSymmetry::Momentum;
		}
		return BrokenSymmetry(sym);
	}

	void SquareGeneral::fill_block_M(int i, int j)
	{
		const OrderType order = this->model->get_order();
		const BrokenSymmetry broken_symmetry = breaks_symmetry(i, j);
		if ((broken_symmetry & BrokenSymmetry::ParticleNumber) && !(order & OrderType::SC)) return;
		if ((broken_symmetry & BrokenSymmetry::Momentum) && !(order & (OrderType::CDW | OrderType::AFM))) return; 

		for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
			for (int k = 0; k < Constants::BASIS_SIZE; k++)
			{
				if (term.delta_momenta.size() > 0U) {
					int l{ k };
					if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
						l = addQTo(k);
					}
					if (term.delta_momenta[0].second.momentum_list.size() == 2U) {
						l = add_index_momenta(l, this->mode_momentum, term.delta_momenta[0].second.momentum_list[1].first);
					}
					M(get_index(i, k), get_index(j, l)) += computeTerm(term, k, l);
				}
				else {
					for (int l = 0; l < Constants::BASIS_SIZE; l++)
					{
						M(get_index(i, k), get_index(j, l)) += computeTerm(term, k, l);
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
						l = add_index_momenta(l, this->mode_momentum, term.delta_momenta[0].second.momentum_list[1].first);
					}
					N(get_index(i, k), get_index(j, l)) += computeTerm(term, k, l);
				}
				else {
					throw std::runtime_error("Offdiagonal term in N!");
					//for (int l = 0; l < Constants::BASIS_SIZE; l++)
					//{
					//	N(get_index(i, k), get_index(j, l)) += computeTerm(term, k, l);
					//}
				}
			}
		}
	}
	void SquareGeneral::setNewModelParameters(mrock::Utility::InputFileReader& input, const Models::ModelParameters& modelParameters)
	{
		this->internal_setNewModelParameters(input, modelParameters);
	}
}