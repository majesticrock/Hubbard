#define _USE_MATH_DEFINES

#include "HubbardCDW.hpp"

namespace Hubbard {
	HubbardCDW::HubbardCDW(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at)
		: BasicHubbardModel(_params, _number_of_basis_terms, _start_basis_at), V(_params.V)
	{
		this->delta_cdw = std::abs(U - V) * 0.5;
		this->delta_sc = std::abs(U + V) * 0.5;
		if (V > 0) {
			this->delta_sc *= 0.25;
		}
		else {
			this->delta_cdw *= 0.25;
		}
		this->delta_eta = 0.01;

		this->hamilton = Matrix_L::Zero(4, 4);
	}
}