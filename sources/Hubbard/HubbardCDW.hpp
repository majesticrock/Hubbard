#pragma once
#include "BasicHubbardModel.hpp"

namespace Hubbard {
	typedef std::complex<double_prec> complex_prec;
	class HubbardCDW : public Model
	{
	protected:
		const complex_prec I = { 0, 1 };

		double_prec V;
		typedef Eigen::Matrix<complex_prec, 4, 4> Matrix_4cL;
		Matrix_4cL complex_h;

		virtual void computeChemicalPotential() override;
		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;

		virtual inline void setParameters(double_prec cdw, double_prec sc, double_prec eta, double_prec cos_n) {
			cdw *= (this->U - this->V) / BASIS_SIZE;
			sc *= this->U / BASIS_SIZE;
			eta *= this->U / BASIS_SIZE;
			cos_n *= (V / (8 * BASIS_SIZE));
			this->delta_cdw = 0.5 * (cdw + this->delta_cdw);
			this->delta_sc = 0.5 * (sc + this->delta_sc);
			this->delta_eta = 0.5 * (eta + this->delta_eta);
			this->delta_occupation = 0.5 * (cos_n + this->delta_occupation);
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