#pragma once
#include "BasicHubbardModel.hpp"

namespace Hubbard {
	class HubbardCDW : public Model
	{
	private:
		typedef Eigen::Vector<complex_prec, 16> ParameterVector;
		inline void printAsRow(ParameterVector& printer) const {
			for (size_t i = 0; i < printer.size(); i++)
			{
				std::cout << " \t" << printer(i);
				if ((i + 1) % 4 == 0) {
					std::cout << "\n\t    ";
				}
			}
			std::cout << std::endl;
		}
	protected:
		double_prec V;
		double_prec V_OVER_N;

		complex_prec delta_sc, delta_cdw, delta_afm, delta_eta;
		complex_prec gamma_sc, xi_sc;
		complex_prec gamma_cdw, xi_cdw;
		complex_prec gamma_afm, xi_afm;
		complex_prec gamma_eta, xi_eta;

		complex_prec delta_occupation_up, delta_occupation_down;
		complex_prec delta_occupation_up_y, delta_occupation_down_y;

		inline double_prec gamma_tilde(double_prec k_x, double_prec k_y) {
			return sin(k_x) + sin(k_y);
		};
		inline double_prec xi_tilde(double_prec k_x, double_prec k_y) {
			return sin(k_x) - sin(k_y);
		};
		virtual void computeChemicalPotential() override;
		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const override {
			return -2. * (
				(delta_occupation_up + 1.) * cos(k_x) + (delta_occupation_up_y + 1.) * cos(k_y)
				).real();
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const override {
			return -2. * (
				(delta_occupation_down + 1.) * cos(k_x) + (delta_occupation_down_y + 1.) * cos(k_y)
				).real();
		};

		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;

		virtual inline void setParameters(ParameterVector& F) {
			F(0)  *= 0.5 * U_OVER_N - 4. * V_OVER_N; // CDW
			F(1)  *= 0.5 * U_OVER_N; // AFM
			F(2)  *= U_OVER_N; // SC
			F(3)  *= V_OVER_N; // Gamma SC
			F(4)  *= V_OVER_N; // Xi SC
			F(5)  *= U_OVER_N; // Eta
			F(6)  *= V_OVER_N; // Occupation Up
			F(7)  *= V_OVER_N; // Occupation Down
			F(8)  *= 0.5 * V_OVER_N; // Gamma CDW
			F(9)  *= 0.5 * V_OVER_N; // Xi CDW
			F(10) *= 0.5 * V_OVER_N; // Gamma AFM
			F(11) *= 0.5 * V_OVER_N; // Xi AFM
			F(12) *= V_OVER_N; // Occupation Up y
			F(13) *= V_OVER_N; // Occupation Down y
			F(14) *= V_OVER_N; // Gamma SC
			F(15) *= V_OVER_N; // Xi SC

			constexpr double_prec new_weight = 0.5;

			this->delta_cdw				  = new_weight * F(0)  + (1 - new_weight) * this->delta_cdw;
			this->delta_afm				  = new_weight * F(1)  + (1 - new_weight) * this->delta_afm;
			this->delta_sc				  = new_weight * F(2)  + (1 - new_weight) * this->delta_sc;
			this->gamma_sc				  = new_weight * F(3)  + (1 - new_weight) * this->gamma_sc;
			this->xi_sc					  = new_weight * F(4)  + (1 - new_weight) * this->xi_sc;
			this->delta_eta				  = new_weight * F(5)  + (1 - new_weight) * this->delta_eta;
			this->delta_occupation_up	  = new_weight * F(6)  + (1 - new_weight) * this->delta_occupation_up;
			this->delta_occupation_down	  = new_weight * F(7)  + (1 - new_weight) * this->delta_occupation_down;
			this->gamma_cdw				  = new_weight * F(8)  + (1 - new_weight) * this->gamma_cdw;
			this->xi_cdw				  = new_weight * F(9)  + (1 - new_weight) * this->xi_cdw;
			this->gamma_afm				  = new_weight * F(10) + (1 - new_weight) * this->gamma_afm;
			this->xi_afm				  = new_weight * F(11) + (1 - new_weight) * this->xi_afm;
			this->delta_occupation_up_y	  = new_weight * F(12) + (1 - new_weight) * this->delta_occupation_up_y;
			this->delta_occupation_down_y = new_weight * F(13) + (1 - new_weight) * this->delta_occupation_down_y;
			this->gamma_eta				  = new_weight * F(14) + (1 - new_weight) * this->gamma_eta;
			this->xi_eta				  = new_weight * F(15) + (1 - new_weight) * this->xi_eta;
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
		HubbardCDW(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at);
		Model::data_set computePhases(const bool print = false) override;
	};
}