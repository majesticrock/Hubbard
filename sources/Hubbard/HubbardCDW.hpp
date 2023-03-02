#pragma once
#include "BasicHubbardModel.hpp"

namespace Hubbard {
	class HubbardCDW : public BasicHubbardModel
	{
	protected:
		double V;

		virtual inline void setParameters(double cdw, double sc, double eta) override {
			this->delta_cdw = ((this->U - this->V) / BASIS_SIZE) * cdw;
			this->delta_sc = ((this->U) / BASIS_SIZE) * sc;
			this->delta_eta = ((this->U) / BASIS_SIZE) * eta;
		};

	public:
		HubbardCDW(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at);
	};
}