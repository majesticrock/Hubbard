#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class UsingBroyden : public Model
	{
	private:
		typedef Eigen::Vector<double_prec, 8> ParameterVector;
	protected:
		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;

		virtual inline void setParameters(ParameterVector& F) {
			constexpr double_prec new_weight = 0.5;

			this->delta_cdw = new_weight * F(0) + (1 - new_weight) * this->delta_cdw;
			this->delta_afm = new_weight * F(1) + (1 - new_weight) * this->delta_afm;
			this->delta_sc = new_weight * F(2) + (1 - new_weight) * this->delta_sc;
			this->gamma_sc = new_weight * F(3) + (1 - new_weight) * this->gamma_sc;
			this->xi_sc = new_weight * F(4) + (1 - new_weight) * this->xi_sc;
			this->delta_eta = new_weight * F(5) + (1 - new_weight) * this->delta_eta;
			this->delta_occupation_up = new_weight * F(6) + (1 - new_weight) * this->delta_occupation_up;
			this->delta_occupation_down = new_weight * F(7) + (1 - new_weight) * this->delta_occupation_down;
		};
	public:
		UsingBroyden(const ModelParameters& _params);

		data_set computePhases(const bool print = false) override;
	};
}