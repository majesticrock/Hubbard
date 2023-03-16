#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class UsingBroyden : public Model
	{
	private:
		Eigen::Vector4d new_sc;
		Eigen::Vector4d new_eta;
		Eigen::Vector4d old_sc;
		Eigen::Vector4d old_eta;
	protected:
		double V;

		virtual void fillHamiltonian(double k_x, double k_y) override;
		virtual inline void setParameters(Eigen::VectorXd& F) {
			F(0) *= (this->U - 8 * this->V) / BASIS_SIZE;
			F(1) *= (this->U) / BASIS_SIZE; // + this->V
			F(2) *= (this->U) / BASIS_SIZE; // + this->V
			this->delta_cdw = 0.5 * (F(0) + this->delta_cdw);
			this->delta_sc = 0.5 * (F(1) + this->delta_sc);
			this->delta_eta = 0.5 * (F(2) + this->delta_eta);
		};
		virtual inline double computeCoefficient(const SymbolicOperators::Coefficient& coeff, const std::optional<Eigen::Vector2i>& momentum = std::nullopt) const override {	
			if (coeff.name == "\\tilde{V}") {
				if (!(momentum.has_value())) throw std::invalid_argument("Calling V without specifying a momentum!");
				// Eventuell ein Faktor 2?
				return (V / BASIS_SIZE) * ( cos( momentum.value()(0) ) + cos( momentum.value()(1) ) );
			}

			return Model::computeCoefficient(coeff, momentum);
		};
		virtual void fill_M_N() override;
	public:
		UsingBroyden(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at);

		data_set computePhases(const bool print = false) override;
	};
}