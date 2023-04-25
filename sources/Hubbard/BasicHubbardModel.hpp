#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class BasicHubbardModel : public Model
	{
	protected:
		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;
		virtual inline void setParameters(double_prec cdw, double_prec sc, double_prec eta) {
			this->delta_cdw = cdw * this->U / BASIS_SIZE;
			this->delta_sc = sc * this->U / BASIS_SIZE;
			this->delta_eta = eta * this->U / BASIS_SIZE;
		};
	public:
		BasicHubbardModel(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at);

		data_set computePhases(const bool print = false) override;
	};
}