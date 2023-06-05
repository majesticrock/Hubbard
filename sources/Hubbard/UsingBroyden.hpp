#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class UsingBroyden : public Model
	{
	private:
		typedef Eigen::Vector<double_prec, 8> ParameterVector;
		inline void printAsRow(const ParameterVector& printer) const {
			for (size_t i = 0; i < printer.size(); i++)
			{
				std::cout << "\t" << printer(i);
				if ((i + 1) % 8 == 0) {
					std::cout << "\n\t    ";
				}
			}
			std::cout << std::endl;
		}

	protected:
		double_prec V;
		double_prec V_OVER_N;

		virtual void computeChemicalPotential() override;
		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const override {
			return -2 * (1. + delta_occupation_up) * gamma(k_x, k_y);
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const override {
			return -2 * (1. + delta_occupation_down) * gamma(k_x, k_y);
		};

		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;

		virtual inline void setParameters(ParameterVector& F) {
			F(0) *= 0.5 * U_OVER_N - 4 * V_OVER_N; // CDW
			F(1) *= 0.5 * U_OVER_N; // AFM
			F(2) *= U_OVER_N; // SC
			F(3) *= V_OVER_N; // Gamma SC
			F(4) *= V_OVER_N; // Xi SC y
			F(5) *= U_OVER_N; // Eta
			F(6) *= V_OVER_N; // Occupation Up
			F(7) *= V_OVER_N; // Occupation Down

			constexpr double_prec new_weight = 0.5;

			this->delta_cdw			    = new_weight * F(0) + (1 - new_weight) * this->delta_cdw;
			this->delta_afm			    = new_weight * F(1) + (1 - new_weight) * this->delta_afm;
			this->delta_sc			    = new_weight * F(2) + (1 - new_weight) * this->delta_sc;
			this->gamma_sc			    = new_weight * F(3) + (1 - new_weight) * this->gamma_sc;
			this->xi_sc				    = new_weight * F(4) + (1 - new_weight) * this->xi_sc;
			this->delta_eta			    = new_weight * F(5) + (1 - new_weight) * this->delta_eta;
			this->delta_occupation_up   = new_weight * F(6) + (1 - new_weight) * this->delta_occupation_up;
			this->delta_occupation_down = new_weight * F(7) + (1 - new_weight) * this->delta_occupation_down;
		};
		virtual inline double_prec computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector2i& momentum) const override {
			if (coeff.name == "\\tilde{V}") {
				//if (!(momentum.has_value())) throw std::invalid_argument("Calling V without specifying a momentum!");
				// Eventuell ein Faktor 2?
				return V_OVER_N * (cos(index_to_k_vector(momentum(0))) + cos(index_to_k_vector(momentum(1))));
			}

			return Model::computeCoefficient(coeff, momentum);
		};
	public:
		UsingBroyden(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at);

		data_set computePhases(const bool print = false) override;
	};
}