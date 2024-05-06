#pragma once
#include "GeneralBasis.hpp"
#include "TermWithDOS.hpp"

namespace Hubbard::Helper {
	template <class DOS>
	class DOSGeneral : public TermWithDOS<DOS>, public GeneralBasis
	{
	private:
		virtual void fill_block_M(int i, int j) override
		{
			// fill M
			for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
				for (int k = 0; k < Constants::BASIS_SIZE; ++k)
				{
					if (!term.delta_momenta.empty()) {
						int l{ k };
						if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
							// delta gamma, -gamma'
							l = this->model->shiftByQ(k);
						}
						M(i + k * number_of_basis_terms, j + l * number_of_basis_terms)
							+= this->computeTerm(term, k, l) * this->approximate_dos[k];
					}
					else {
						for (int l = 0; l < Constants::BASIS_SIZE; ++l)
						{
							M(i + k * number_of_basis_terms, j + l * number_of_basis_terms)
								+= this->computeTerm(term, k, l) * (this->approximate_dos[k] * this->approximate_dos[l] / this->INV_GAMMA_DISC);
						}
					}
				}
			}
		};

		virtual void fill_block_N(int i, int j) override
		{
			// fill N
			for (const auto& term : wicks.N[number_of_basis_terms * j + i]) {
				for (int k = 0; k < Constants::BASIS_SIZE; ++k)
				{
					if (!term.delta_momenta.empty()) {
						int l{ k };
						if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
							// delta gamma, -gamma'
							l = this->model->shiftByQ(k);
						}
						N(i + k * number_of_basis_terms, j + l * number_of_basis_terms)
							+= this->computeTerm(term, k, l) * this->approximate_dos[k];
					}
					else {
						for (int l = 0; l < Constants::BASIS_SIZE; ++l)
						{
							N(i + k * number_of_basis_terms, j + l * number_of_basis_terms)
								+= this->computeTerm(term, k, l) * (this->approximate_dos[k] * this->approximate_dos[l] / this->INV_GAMMA_DISC);
						}
					}
				}
			}
		}

	public:
		virtual const BaseModel<global_floating_type>& getModel() const override {
			return *(this->model);
		};
		virtual BaseModel<global_floating_type>& getModel() override {
			return *(this->model);
		};

		DOSGeneral(Utility::InputFileReader& input, const ModelParameters& modelParameters) : TermWithDOS<DOS>(input, modelParameters), GeneralBasis(input) {
			this->TOTAL_BASIS = this->number_of_basis_terms * Constants::BASIS_SIZE;
		};

		void setNewModelParameters(Utility::InputFileReader& input, const ModelParameters& modelParameters) override {
			this->internal_setNewModelParameters(input, modelParameters);
		};
	};
}