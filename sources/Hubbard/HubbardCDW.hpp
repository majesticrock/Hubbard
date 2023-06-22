#pragma once
#include "Model.hpp"

namespace Hubbard {
	class HubbardCDW : public Model
	{
	private:
		typedef Eigen::Vector<complex_prec, 16> ParameterVector;
	protected:
		std::vector<complex_prec*> param_mapper;
		std::vector<double_prec> param_coefficients;

		complex_prec c_delta_sc, c_delta_cdw, c_delta_afm, c_delta_eta;
		complex_prec c_gamma_sc, c_xi_sc;
		complex_prec c_gamma_cdw, c_xi_cdw;
		complex_prec c_gamma_afm, c_xi_afm;
		complex_prec c_gamma_eta, c_xi_eta;

		complex_prec c_delta_occupation_up, c_delta_occupation_down;
		complex_prec c_delta_occupation_up_y, c_delta_occupation_down_y;

		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const override {
			return -2. * (
				(c_delta_occupation_up + 1.) * cos(k_x) + (c_delta_occupation_up_y + 1.) * cos(k_y)
				).real();
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const override {
			return -2. * (
				(c_delta_occupation_down + 1.) * cos(k_x) + (c_delta_occupation_down_y + 1.) * cos(k_y)
				).real();
		};

		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;

		void addToParameterSet(const SpinorMatrix& rho, Eigen::Ref<ParameterVector> F, double_prec k_x, double_prec k_y) {
			F(0) -= (rho(0, 1) + rho(1, 0) - rho(2, 3) - rho(3, 2)).real(); // CDW
			F(1) -= (rho(0, 1) + rho(1, 0) + rho(2, 3) + rho(3, 2)).real(); // AFM
			F(2) -= (rho(0, 2) + rho(1, 3)); // SC
			F(3) -= gamma(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Gamma SC
			F(4) -= xi(k_x, k_y) * (rho(0, 2) - rho(1, 3)); // Xi SC
			F(5) -= (rho(0, 3) + rho(1, 2)); // Eta
			F(6) -= cos(k_x) * (rho(0, 0) - rho(1, 1)).real(); // Occupation Up
			F(7) += cos(k_x) * (rho(2, 2) - rho(3, 3)).real(); // Occupation Down
			F(8) += I * gamma(k_x, k_y) * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Gamma CDW
			F(9) += I * xi(k_x, k_y) * (rho(0, 1) - rho(1, 0) + rho(2, 3) - rho(3, 2)).imag(); // Xi CDW
			F(10) += I * gamma(k_x, k_y) * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Gamma AFM
			F(11) += I * xi(k_x, k_y) * (rho(0, 1) - rho(1, 0) - rho(2, 3) + rho(3, 2)).imag(); // Xi AFM
			F(12) -= cos(k_y) * (rho(0, 0) - rho(1, 1)).real(); // Occupation Up y
			F(13) += cos(k_y) * (rho(2, 2) - rho(3, 3)).real(); // Occupation Down y
			F(14) -= gamma(k_x, k_y) * (rho(0, 3) - rho(1, 2)); // Gamma eta
			F(15) -= xi(k_x, k_y) * (rho(0, 3) - rho(1, 2)); // Xi eta
		};

		template <const int vector_size>
		void setParameters(Eigen::Vector<complex_prec, vector_size>& F) {
			for (size_t i = 0; i < vector_size; i++)
			{
				F(i) *= param_coefficients[i];
			}

			constexpr double_prec new_weight = 0.5;
			for (size_t i = 0; i < vector_size; i++)
			{
				*(param_mapper[i]) = new_weight * F(i) + (1 - new_weight) * (*(param_mapper[i]));
			}
		};
	public:
		HubbardCDW(const ModelParameters& _params);
		virtual Model::data_set computePhases(const bool print = false) override;
	};
}