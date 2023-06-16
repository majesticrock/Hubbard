#pragma once
#include "Model.hpp"

namespace Hubbard {
	constexpr size_t NUMBER_OF_PARAMETERS = 16;

	class HubbardCDW : public Model
	{
	private:
		typedef Eigen::Vector<complex_prec, NUMBER_OF_PARAMETERS> ParameterVector;
		inline void printAsRow(const ParameterVector& printer) const {
			for (size_t i = 0; i < printer.size(); i++)
			{
				std::cout << " \t" << printer(i);
				if ((i + 1) % 4 == 0) {
					std::cout << "\n\t    ";
				}
			}
			std::cout << std::endl;
		}

		std::vector<complex_prec*> param_mapper;
		std::vector<double_prec> param_coefficients;
	protected:
		complex_prec c_delta_sc,  c_delta_cdw, c_delta_afm, c_delta_eta;
		complex_prec c_gamma_sc,  c_xi_sc;
		complex_prec c_gamma_cdw, c_xi_cdw;
		complex_prec c_gamma_afm, c_xi_afm;
		complex_prec c_gamma_eta, c_xi_eta;

		complex_prec c_delta_occupation_up, c_delta_occupation_down;
		complex_prec c_delta_occupation_up_y, c_delta_occupation_down_y;

		inline double_prec gamma_tilde(double_prec k_x, double_prec k_y) {
			return sin(k_x) + sin(k_y);
		};
		inline double_prec xi_tilde(double_prec k_x, double_prec k_y) {
			return sin(k_x) - sin(k_y);
		};
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

		virtual inline void setParameters(ParameterVector& F) {
			for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
			{
				F(i) *= param_coefficients[i];
			}

			constexpr double_prec new_weight = 0.5;
			for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
			{
				*(param_mapper[i]) = new_weight * F(i) + (1 - new_weight) * (*(param_mapper[i]));
			}
		};

	public:
		HubbardCDW(const ModelParameters& _params);
		Model::data_set computePhases(const bool print = false) override;
	};
}