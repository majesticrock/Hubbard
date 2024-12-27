#pragma once
#include "XPModes.hpp"
#include "TermWithDOS.hpp"

#ifdef _EXACT_DOS
#define _FACTOR_DIAG (this->approximate_dos[k] * this->getWeightFromIndex(k))
#define _FACTOR_OFF_DIAG (this->approximate_dos[k] * this->approximate_dos[l] * this->getWeightFromIndex(k) * this->getWeightFromIndex(l))
#else
#define _FACTOR_DIAG this->approximate_dos[k]
#define _FACTOR_OFF_DIAG (this->approximate_dos[k] * this->approximate_dos[l] / this->INV_GAMMA_DISC)
#endif

namespace Hubbard::Helper {
	template <class DOS>
	class DOS_XP : public TermWithDOS<DOS>, public XPModes
	{
	protected:
		inline global_floating_type computeRealTerm(const mrock::SymbolicOperators::WickTerm& term, int k, int l) const {
			const auto result = this->computeTerm(term, k, l);
			if (abs(std::imag(result)) > ERROR_MARGIN) {
				throw std::runtime_error("computeRealTerm() encountered a complex value!");
			}
			return result;
		};

		// K_+ / K_-
		virtual void fill_block_M(int i, int j) override
		{
			const int sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()
				? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;
			const int inner_sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), j) == cdw_basis_positions.end()
				? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;

			for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
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

						if (i < hermitian_size) {
							K_plus(hermitian_offsets[i] + k, hermitian_offsets[j] + l)
								+= computeRealTerm(term, k, l) * _FACTOR_DIAG;
						}
						else {
							K_minus(antihermitian_offsets[i - hermitian_size] + k, antihermitian_offsets[j - hermitian_size] + l)
								+= computeRealTerm(term, k, l) * _FACTOR_DIAG;
						}
					}
					else {
						for (int l = 0; l < inner_sum_limit; ++l)
						{
							if (i < hermitian_size) {
								K_plus(hermitian_offsets[i] + k, hermitian_offsets[j] + l)
									+= computeRealTerm(term, k, l) * _FACTOR_OFF_DIAG;
							}
							else {
								K_minus(antihermitian_offsets[i - hermitian_size] + k, antihermitian_offsets[j - hermitian_size] + l)
									+= computeRealTerm(term, k, l) * _FACTOR_OFF_DIAG;
							}
						}
					}
				} // end k-loop
			} // end term-loop
		}
		// L
		virtual void fill_block_N(int i, int j) override
		{
			const int sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), i) == cdw_basis_positions.end()
				? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;
			const int inner_sum_limit = std::find(cdw_basis_positions.begin(), cdw_basis_positions.end(), j) == cdw_basis_positions.end()
				? Constants::BASIS_SIZE : Constants::BASIS_SIZE / 2;

			for (const auto& term : wicks.N[number_of_basis_terms * j + i]) {
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

						L(hermitian_offsets[i] + k, antihermitian_offsets[j - hermitian_size] + l)
							+= computeRealTerm(term, k, l) * _FACTOR_DIAG;
						// Technically, we need to multiply this term by h(gamma - gamma') = Delta gamma
						// But we dont, that's why we need to devide it later in the offdiagonal parts
						// The factor is added, when we compute the lanczos coefficients
						// If we were to include it here, the matrix elements would become very small (as Delta gamma is small)
						// which is bad for numerical stability
					}
					else {
						for (int l = 0; l < inner_sum_limit; l++)
						{
							L(hermitian_offsets[i] + k, antihermitian_offsets[j - hermitian_size] + l)
								+= computeRealTerm(term, k, l) * _FACTOR_OFF_DIAG;
						}
					}
				} // end k-loop
			} // end term-loop
		}

	public:
		virtual const Models::BaseModel<global_floating_type>& getModel() const override {
			return *(this->model);
		};
		virtual Models::BaseModel<global_floating_type>& getModel() override {
			return *(this->model);
		};

		DOS_XP(mrock::Utility::InputFileReader& input, const Models::ModelParameters& modelParameters) : TermWithDOS<DOS>(input, modelParameters), XPModes(input) {
			this->TOTAL_BASIS = this->number_of_basis_terms * Constants::BASIS_SIZE;
		};

		void setNewModelParameters(mrock::Utility::InputFileReader& input, const Models::ModelParameters& modelParameters) override {
			this->internal_setNewModelParameters(input, modelParameters);
		};
	};
}

#undef _FACTOR_DIAG
#undef _FACTOR_OFF_DIAG