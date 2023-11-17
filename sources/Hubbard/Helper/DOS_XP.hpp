#pragma once
#include "XPModes.hpp"
#include "TermWithDOS.hpp"

namespace Hubbard::Helper {
	template <class DOS>
	class DOS_XP : public XPModes, public TermWithDOS<DOS>
	{
	private:
		inline global_floating_type computeRealTerm(const SymbolicOperators::WickTerm& term, int k, int l) const {
			const auto result = this->computeTerm(term, k, l);
			if (abs(result.imag()) > ERROR_MARGIN) {
				throw std::runtime_error("computeRealTerm() encountered a complex value!");
			}
			return result.real();
		};

		virtual void fillBlock(int i, int j) override
		{
			constexpr std::array<int, 4> cdw_basis_positions{ 2,3,8,9 };
			const int hermitian_offsets[6] = {
				0,							Constants::BASIS_SIZE,
				2 * Constants::BASIS_SIZE,	(5 * Constants::BASIS_SIZE) / 2,
				3 * Constants::BASIS_SIZE,	4 * Constants::BASIS_SIZE
			};
			const int antihermitian_offsets[4] = {
				0,									Constants::BASIS_SIZE,
				2 * Constants::BASIS_SIZE,			(5 * Constants::BASIS_SIZE) / 2
			};

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

	public:
		virtual const BaseModel<global_floating_type>& getModel() const override {
			return *(this->model);
		};
		virtual BaseModel<global_floating_type>& getModel() override {
			return *(this->model);
		};

		DOS_XP(Utility::InputFileReader& input) : XPModes(input), TermWithDOS<DOS>(input) {
			this->TOTAL_BASIS = this->number_of_basis_terms * Constants::BASIS_SIZE;
		};
	};
}