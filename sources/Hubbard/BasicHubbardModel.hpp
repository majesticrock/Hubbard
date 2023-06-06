#pragma once
#include "Model.hpp"
#include "Constants.hpp"

namespace Hubbard {
	class BasicHubbardModel : public Model
	{
	private:
		typedef Eigen::Vector<double_prec, 4> ParameterVector;
		inline void printAsRow(ParameterVector& printer) const {
			for (size_t i = 0; i < printer.size(); i++)
			{
				std::cout << "\t" << printer(i);
			}
			std::cout << std::endl;
		};
	protected:
		virtual void fillHamiltonian(double_prec k_x, double_prec k_y) override;
		virtual inline void setParameters(ParameterVector& F) {
			F(0) *= 0.5 * this->U / BASIS_SIZE; // CDW
			F(1) *= 0.5 * this->U / BASIS_SIZE; // AFM
			F(2) *= this->U / BASIS_SIZE; // SC
			F(3) *= this->U / BASIS_SIZE; // Eta

			this->delta_cdw = 0.5 * (F(0) + this->delta_cdw);
			this->delta_afm = 0.5 * (F(1) + this->delta_afm);
			this->delta_sc = 0.5 * (F(2) + this->delta_sc);
			this->delta_eta = 0.5 * (F(3) + this->delta_eta);
		};
	public:
		BasicHubbardModel(const ModelParameters& _params, int _number_of_basis_terms, int _start_basis_at);

		data_set computePhases(const bool print = false) override;
	};
}