#pragma once
#include "DetailModelConstructor.hpp"
#include "../DOSModels/BroydenDOS.hpp"

namespace Hubbard::Helper {
	using DOS = DensityOfStates::Square;

	class TermWithDOS : protected DetailModelConstructor<DOSModels::BroydenDOS<DOS>>
	{
	private:
		using Integrator = DOS::Integrator<complex_prec>;
		Integrator _integrator{};

	protected:
		complex_prec getExpectationValue(const SymbolicOperators::WickOperator& op, int gamma_idx) const {
			auto it = wick_map.find(op.type);
			if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + op.type);

			int index = it->second;
			if (op.type == "g" || op.type == "n") {
				auto jt = wick_spin_offset.find(op.indizes[0]);
				if (jt == wick_map.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
				index += jt->second;
			}

			if (op.isDaggered) return std::conj(expecs[index](gamma_idx, 0));
			return expecs[index](gamma_idx, 0);
		};

		complex_prec computeTerm(const SymbolicOperators::WickTerm& term, int gamma_idx, int gamma_prime_idx) const {
			checkTermValidity(term);
			const global_floating_type gamma{ DOS::abscissa_v(gamma_idx) };
			const global_floating_type gamma_prime{ DOS::abscissa_v(gamma_prime_idx) };

			auto getCoefficient = [&](global_floating_type x) {
				return this->model->computeCoefficient(term.coefficients.back(), x);
			};

			auto getCoefficientAndExpec = [&](global_floating_type x, size_t expec_pos) {
				return this->model->computeCoefficient(term.coefficients.back(), x) * getExpectationValue(term.operators[expec_pos], x);
			};

			if (term.isIdentity()) {
				if (term.hasSingleCoefficient()) {
					return term.multiplicity * getCoefficient(gamma);
				}
				return static_cast<global_floating_type>(term.multiplicity);
			}
			if (term.sum_momenta.size() > 0U) {
				if (term.isBilinear()) {
					// bilinear
					if (term.hasSingleCoefficient()) {
						if (term.coefficients.back().momentum.isUsed('q') != -1) {
							if (term.operators.front().momentum.isUsed('q') == -1) return 0;

							auto _integrate_function = [&](size_t index) -> complex_prec {
								return DOS::abscissa_v(index) *	
									(getCoefficientAndExpec(DOS::abscissa_v(index), 0U) - getCoefficientAndExpec(DOS::abscissa_v(2 * index), 0U));
							};
							return term.multiplicity * getCoefficient(gamma) * _integrator.integrate_by_index(_integrate_function);
						}
						else {
							return term.multiplicity * getCoefficient(gamma) * getSumOfAll(term.operators[0]);
						}
					}
					else {
					}
				}
				else if (term.isQuartic()) {
					// quartic
				}
			}

			throw std::runtime_error("Not yet implemented!");
		};

	public:
		TermWithDOS(Utility::InputFileReader& input) : DetailModelConstructor(input) {};
	};
}