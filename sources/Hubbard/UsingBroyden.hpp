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
		double cos_occupation;
		double cos_occupation_old;
	protected:
		long double V;

		virtual void computeChemicalPotential() override;
		virtual void fillHamiltonian(double k_x, double k_y) override;
		virtual inline void setParameters(Eigen::VectorXd& F) {
			F(0) *= (this->U - this->V) / BASIS_SIZE;
			F(1) *= (this->U) / BASIS_SIZE; // + this->V
			F(2) *= (this->U) / BASIS_SIZE; // + this->V
			this->delta_cdw = 0.5 * (F(0) + this->delta_cdw);
			this->delta_sc = 0.5 * (F(1) + this->delta_sc);
			this->delta_eta = 0.5 * (F(2) + this->delta_eta);

			this->cos_occupation *= ((-V) / (16 * BASIS_SIZE));
			this->cos_occupation_old = cos_occupation;
		};
		virtual inline long double computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector2i& momentum) const override {
			if (coeff.name == "\\tilde{V}") {
				//if (!(momentum.has_value())) throw std::invalid_argument("Calling V without specifying a momentum!");
				// Eventuell ein Faktor 2?
				return (V / (16 * BASIS_SIZE)) * 2 * (cos(index_to_k_vector(momentum(0))) + cos(index_to_k_vector(momentum(1))));
			}

			return Model::computeCoefficient(coeff, momentum);
		};
	public:
		UsingBroyden(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at);

		data_set computePhases(const bool print = false) override;
	};
}