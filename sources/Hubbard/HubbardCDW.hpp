#pragma once
#include "BasicHubbardModel.hpp"

namespace Hubbard {
	typedef std::complex<double_prec> complex_prec;
	class HubbardCDW : public Model
	{
	protected:
		typedef Eigen::Vector<double_prec, 10> ParameterVector;
		typedef Eigen::Matrix<complex_prec, 4, 4> Matrix_4cL;
		const complex_prec I = { 0, 1 };

		double_prec V;
		double_prec xi_sc_x, xi_sc_y, xi_eta_x, xi_eta_y;
		Matrix_4cL complex_h;

		virtual void computeChemicalPotential() override;
		inline virtual double_prec renormalizedEnergy_up(double_prec k_x, double_prec k_y) const override {
			return -2 * (1. + delta_occupation_up) * (cos(k_x) + cos(k_y));
		};
		inline virtual double_prec renormalizedEnergy_down(double_prec k_x, double_prec k_y) const override {
			return -2 * (1. + delta_occupation_down) * (cos(k_x) + cos(k_y));
		};

		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;

		virtual inline void setParameters(ParameterVector& F) {
			//F(0) *= (this->U - this->V) / BASIS_SIZE; // CDW Up
			//F(1) *= (this->U - this->V) / BASIS_SIZE; // CDW Down
			auto buf_up = ((this->U - 0.5 * this->V) * F(0) - 0.5 * this->V * F(1)) / BASIS_SIZE;
			auto buf_down = ((this->U - 0.5 * this->V) * F(1) - 0.5 * this->V * F(0)) / BASIS_SIZE;
			F(0) = buf_up;
			F(1) = buf_down;
			F(2) *= this->U / BASIS_SIZE; // SC
			F(3) *= this->U / BASIS_SIZE; // Eta
			F(4) *= V / (8 * BASIS_SIZE); // Occupation Up
			F(5) *= V / (8 * BASIS_SIZE); // Occupation Down
			F(6) *= V / (8 * BASIS_SIZE); // Xi SC x
			F(7) *= V / (8 * BASIS_SIZE); // Xi SC y
			F(8) *= V / (8 * BASIS_SIZE); // Xi Eta x
			F(9) *= V / (8 * BASIS_SIZE); // Xi Eta y

			this->delta_cdw_up = 0.5 * (F(0) + this->delta_cdw_up);
			this->delta_cdw_down = 0.5 * (F(1) + this->delta_cdw_down);
			this->delta_sc = 0.5 * (F(2) + this->delta_sc);
			this->delta_eta = 0.5 * (F(3) + this->delta_eta);
			this->delta_occupation_up = 0.5 * (F(4) + this->delta_occupation_up);
			this->delta_occupation_down = 0.5 * (F(5) + this->delta_occupation_down);
			this->xi_sc_x = 0.5 * (F(6) + this->xi_sc_x);
			this->xi_sc_y = 0.5 * (F(7) + this->xi_sc_y);
			this->xi_eta_x = 0.5 * (F(8) + this->xi_eta_x);
			this->xi_eta_y = 0.5 * (F(9) + this->xi_eta_y);
		};
		virtual inline double_prec computeCoefficient(const SymbolicOperators::Coefficient& coeff, const Eigen::Vector2i& momentum) const override {
			if (coeff.name == "\\tilde{V}") {
				//if (!(momentum.has_value())) throw std::invalid_argument("Calling V without specifying a momentum!");
				// Eventuell ein Faktor 2?
				return (V / (8 * BASIS_SIZE)) * (cos(index_to_k_vector(momentum(0))) + cos(index_to_k_vector(momentum(1))));
			}

			return Model::computeCoefficient(coeff, momentum);
		};

	public:
		HubbardCDW(ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at);
		Model::data_set computePhases(const bool print = false) override;
	};
}