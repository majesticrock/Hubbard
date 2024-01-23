#pragma once
#include "XPModes.hpp"
#include "TermWithDOS.hpp"

namespace Hubbard::Helper {
	template <class DOS>
	class DOS_XP : public TermWithDOS<DOS>, public XPModes
	{
	private:
		inline global_floating_type computeRealTerm(const SymbolicOperators::WickTerm& term, int k, int l) const {
			const auto result = this->computeTerm(term, k, l);
			if (abs(result.imag()) > ERROR_MARGIN) {
				throw std::runtime_error("computeRealTerm() encountered a complex value!");
			}
			return result.real();
		};
		
		virtual void fill_block_M(int i, int j) override
		{
			const int sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()
				? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;
			const int inner_sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), j) == cdw_basis_positions.end()
				? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;

			// K_+ / K_-
			// Ignore the offdiagonal blocks as they are 0
			if (i < 6 && j > 5) return;
			if (j < 6 && i > 5) return;

			for (const auto& term : wicks_M[number_of_basis_terms * j + i]) {
				for (int k = 0; k < sum_limit; ++k)
				{
					if (term.delta_momenta.size() > 0U) {
						int l{ k };
						if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
							l = this->model->shiftByQ(k);
						}
						if (l >= inner_sum_limit) {
							continue;
						}

						if (i < 6) {
							K_plus(hermitian_offsets[i] + k, hermitian_offsets[j] + l)
								+= computeRealTerm(term, k, l) * this->approximate_dos[k];
						}
						else {
							K_minus(antihermitian_offsets[i - 6] + k, antihermitian_offsets[j - 6] + l)
								+= computeRealTerm(term, k, l) * this->approximate_dos[k];
						}
					}
					else {
						for (int l = 0; l < inner_sum_limit; ++l)
						{
							if (i < 6) {
								K_plus(hermitian_offsets[i] + k, hermitian_offsets[j] + l)
									+= computeRealTerm(term, k, l) * (this->approximate_dos[k] * this->approximate_dos[l] / this->INV_GAMMA_DISC);
							}
							else {
								K_minus(antihermitian_offsets[i - 6] + k, antihermitian_offsets[j - 6] + l)
									+= computeRealTerm(term, k, l) * (this->approximate_dos[k] * this->approximate_dos[l] / this->INV_GAMMA_DISC);
							}
						}
					}
				} // end k-loop
			} // end term-loop
		}

		virtual void fill_block_N(int i, int j) override
		{
			const int sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()
				? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;
			const int inner_sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), j) == cdw_basis_positions.end()
				? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;

			// L
			if (i < 6 && j > 5) {
				for (const auto& term : wicks_N[number_of_basis_terms * j + i]) {
					for (int k = 0; k < sum_limit; ++k)
					{
						if (term.delta_momenta.size() > 0U) {
							int l{ k };
							if (term.delta_momenta[0].first.add_Q != term.delta_momenta[0].second.add_Q) {
								l = this->model->shiftByQ(k);
							}
							if (l >= inner_sum_limit) {
								continue;
							}

							L(hermitian_offsets[i] + k, antihermitian_offsets[j - 6] + l)
								+= computeRealTerm(term, k, l) * this->approximate_dos[k];
							// Technically, we need to multiply this term by h(gamma - gamma') = Delta gamma
							// But we dont, that's why we need to devide it later in the offdiagonal parts
							// The factor is added, when we compute the lanczos coefficients
							// If we were to include it here, the matrix elements would become very small (as Delta gamma is small)
							// which is bad for numerical stability
						}
						else {
							for (int l = 0; l < inner_sum_limit; l++)
							{
								L(hermitian_offsets[i] + k, antihermitian_offsets[j - 6] + l)
									+= computeRealTerm(term, k, l) * (this->approximate_dos[k] * this->approximate_dos[l] / this->INV_GAMMA_DISC);
							}
						}
					} // end k-loop
				} // end term-loop
			}
		}

	public:
		virtual const BaseModel<global_floating_type>& getModel() const override {
			return *(this->model);
		};
		virtual BaseModel<global_floating_type>& getModel() override {
			return *(this->model);
		};

		DOS_XP(Utility::InputFileReader& input, const ModelParameters& modelParameters) : TermWithDOS<DOS>(input, modelParameters), XPModes(input) {
			this->TOTAL_BASIS = this->number_of_basis_terms * Constants::BASIS_SIZE;
		};
	};
}