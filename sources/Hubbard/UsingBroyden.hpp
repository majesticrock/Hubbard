#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class UsingBroyden : public Model
	{
	private:
		typedef Eigen::Vector<double_prec, 8> ParameterVector;
		inline void printAsRow(ParameterVector& printer) const {
			for (size_t i = 0; i < printer.size(); i++)
			{
				std::cout << "\t" << printer(i);
				if ((i + 1) % 7 == 0) {
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
			return -2 * (1. + delta_occupation_up) * (cos(k_x) + cos(k_y));
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const override {
			return -2 * (1. + delta_occupation_down) * (cos(k_x) + cos(k_y));
		};

		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;

		virtual inline void setParameters(ParameterVector& F) {
			auto buf_up = U_OVER_N * F(0) - 4 * V_OVER_N * (F(0) + F(1));
			auto buf_down = U_OVER_N * F(1) - 4 * V_OVER_N * (F(0) + F(1));
			F(0) = buf_up;
			F(1) = buf_down;
			F(2) *= U_OVER_N; // SC
			F(3) *= V_OVER_N; // Gamma SC
			F(4) *= V_OVER_N; // Xi SC y
			F(5) *= U_OVER_N; // Eta
			F(6) *= V_OVER_N; // Occupation Up
			F(7) *= V_OVER_N; // Occupation Down

			this->delta_cdw_up = 0.5 * (F(0) + this->delta_cdw_up);
			this->delta_cdw_down = 0.5 * (F(1) + this->delta_cdw_down);
			this->delta_sc = 0.5 * (F(2) + this->delta_sc);
			this->gamma_sc = 0.5 * (F(3) + this->gamma_sc);
			this->xi_sc = 0.5 * (F(4) + this->xi_sc);
			this->delta_eta = 0.5 * (F(5) + this->delta_eta);
			this->delta_occupation_up = 0.5 * (F(6) + this->delta_occupation_up);
			this->delta_occupation_down = 0.5 * (F(7) + this->delta_occupation_down);
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