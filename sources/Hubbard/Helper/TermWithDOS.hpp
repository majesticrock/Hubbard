#pragma once
#include "DetailModelConstructor.hpp"
#include "../DOSModels/BroydenDOS.hpp"
#include "../NumericalMomentum.hpp"
#include <complex>
#include <cmath>

namespace Hubbard::Helper {
	template <class DOS>
	class TermWithDOS : protected DetailModelConstructor<DOSModels::BroydenDOS<DOS>>
	{
	private:
		using Integrator = typename DOS::Integrator<complex_prec>;
		Integrator _integrator{};

		// The non-summing terms (except for the epsilon terms) are proportioal to 1/#lattice sites
		// This number is probably arbitrary in the DOS-description so we will fix it to a value that works
	protected:
		std::vector<global_floating_type> approximate_dos;

		complex_prec getExpectationValue(const SymbolicOperators::WickOperator& op, int gamma_idx) const {
			auto it = this->wick_map.find(op.type);
			if (it == this->wick_map.end()) throw std::invalid_argument("Term type not recognized: " + op.type);

			int index = it->second;
			if (op.type == "g" || op.type == "n") {
				auto jt = this->wick_spin_offset.find(op.indizes[0]);
				if (jt == this->wick_map.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
				index += jt->second;
			}

			auto offset = [&op, &gamma_idx, this]() -> int {
				if (op.momentum.add_Q) {
					return this->model->shiftByQ(gamma_idx);
				}
				return gamma_idx;
				};

			if (op.isDaggered) return std::conj(this->expecs[index](offset(), 0));
			return this->expecs[index](offset(), 0);
		};

		complex_prec computeTerm(const SymbolicOperators::WickTerm& term, int gamma_idx, int gamma_prime_idx) {
			const global_floating_type& gamma{ this->model->getGammaFromIndex(gamma_idx) };
			const global_floating_type& gamma_prime{ this->model->getGammaFromIndex(gamma_prime_idx) };

			auto getCoefficient = [&]() {
				if (term.getFirstCoefficient().dependsOn('l')) {
					return gamma_prime * term.getFactor() * this->model->computeCoefficient(term.getFirstCoefficient(), gamma);
				}
				return term.getFactor() * this->model->computeCoefficient(term.getFirstCoefficient(), gamma);
				};

			auto getCoefficientAndExpec = [&](size_t expec_pos) {
				return getCoefficient() * getExpectationValue(term.operators[expec_pos], gamma_idx);
				};

			if (term.isIdentity()) {
				if (term.hasSingleCoefficient()) {
					return getCoefficient();
				}
				return term.getFactor();
			}
			if (term.sum_momenta.size() > 0U) {
				if (term.isBilinear()) {
					// bilinear
					return getCoefficient() * this->getSumOfAll(term.operators.front(), term.getFirstCoefficient().dependsOn('q'));
				}
				else if (term.isQuartic()) {
					// quartic
					int q_dependend = term.whichOperatorDependsOn('q');
					if (q_dependend < 0) throw std::invalid_argument("Suming q, but no operator depends on q.");

					return getCoefficientAndExpec(q_dependend == 0) * this->getSumOfAll(term.operators[q_dependend], term.getFirstCoefficient().dependsOn('q'));
				}
			}

			if (term.hasSingleCoefficient()) {
				// Can never be an identity (checked above) and only be bilinear or quartic (checked in validity)
				complex_prec ret{};
				if (term.isBilinear()) {
					const auto& op = term.operators.front();
					ret = (op.dependsOn('l') ? getExpectationValue(op, gamma_prime_idx) : getExpectationValue(op, gamma_idx));
				}
				else {
					int l_dependend = term.whichOperatorDependsOn('l');
					ret = getExpectationValue(term.operators[l_dependend == 0], gamma_idx)
						* getExpectationValue(term.operators[l_dependend], gamma_prime_idx);
				}
				return ret * getCoefficient();
			}

			return term.getFactor() * getExpectationValue(term.operators[0U], gamma_idx);
		};

	public:
		TermWithDOS(Utility::InputFileReader& input) : DetailModelConstructor<DOSModels::BroydenDOS<DOS>>(input),
			approximate_dos(Constants::BASIS_SIZE, 0.0)
		{
#ifdef _EXACT_DOS
			for (size_t i = 0U; i < Constants::BASIS_SIZE; ++i)
			{
				approximate_dos[i] = DOS::values_v(i) * DOS::weights_v(i);
			}
#else
			constexpr int faktor = 20;
			Constants::K_DISCRETIZATION *= faktor;
			Constants::PI_DIV_DISCRETIZATION /= faktor;

			NumericalMomentum<DOS::DIMENSION> ks;
			do {
				//std::cout << ks.gamma() << "\t" << 0.5 * (1 - ks.gamma() / DOS::LOWER_BORDER) * (Constants::BASIS_SIZE - 1) << std::endl;
				// Irgendeine Zahl zwischen 0 und BASIS_SIZE
				approximate_dos.at(std::round(0.5 * (1 - ks.gamma() / DOS::LOWER_BORDER) * (Constants::BASIS_SIZE - 1))) += 1;
			} while (ks.iterateFullBZ());
			for (int i = 0U; i < Constants::HALF_BASIS; ++i)
			{
				approximate_dos[i] += approximate_dos[Constants::BASIS_SIZE - 1 - i];
				approximate_dos[i] *= 0.5;
				approximate_dos[Constants::BASIS_SIZE - 1 - i] = approximate_dos[i];
			}
			for (auto& value : this->approximate_dos)
			{
				value *= ((double) Constants::BASIS_SIZE) / (4. * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
			}

			Constants::K_DISCRETIZATION /= faktor;
			Constants::PI_DIV_DISCRETIZATION *= faktor;
#endif
		};
	};
}