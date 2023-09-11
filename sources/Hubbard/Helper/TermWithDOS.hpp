#pragma once
#include "DetailModelConstructor.hpp"
#include "../DOSModels/BroydenDOS.hpp"
#include <complex>

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

		complex_prec computeTerm(const SymbolicOperators::WickTerm& term, int gamma_idx, int gamma_prime_idx) {
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
					return term.getFactor() * getCoefficient(gamma);
				}
				return term.getFactor();
			}
			if (term.sum_momenta.size() > 0U) {
				if (term.isBilinear()) {
					// bilinear
					if (term.coefficients.back().dependsOn('q')) {
						auto _integrate_lambda = [&](size_t index) -> complex_prec {
							return DOS::abscissa_v(index)
								* (getExpectationValue(term.operators[0U], DOS::abscissa_v(index))
									- getExpectationValue(term.operators[0U], DOS::abscissa_v(2 * index)));
						};
						return term.getFactor() * getCoefficient(gamma) * _integrator.integrate_by_index(_integrate_lambda);
					}
					return term.getFactor() * getCoefficient(gamma) * getSumOfAll(term.operators.front());
				}
				else if (term.isQuartic()) {
					// quartic
					int q_dependend = term.whichOperatorDependsOn('q');
					if (q_dependend < 0) throw std::invalid_argument("Suming q, but no operator depends on q.");

					if (term.coefficients.back().dependsOn('q')) {
						auto _integrate_lambda = [&](size_t index) -> complex_prec {
							return DOS::abscissa_v(index)
								* (getExpectationValue(term.operators[q_dependend], DOS::abscissa_v(index))
									- getExpectationValue(term.operators[q_dependend], DOS::abscissa_v(2 * index)));
						};
						// q_dependend can either be 1 or 0; the other operator always depends solely on k
						// Hence q_dependend == 0 gives the positions of the k-dependend operator
						return term.getFactor()
							* getCoefficientAndExpec(gamma, q_dependend == 0) * _integrator.integrate_by_index(_integrate_lambda);
					}
					return term.getFactor()
						* getCoefficientAndExpec(gamma, q_dependend == 0) * getSumOfAll(term.operators.front());
				}
			}

			if (term.hasSingleCoefficient()) {
				if (term.coefficients.back().name != "\\epsilon_0") {
					// This we should be able to choose this as zero?
					// To be verified
					return 0;
				}
				// Can never be an identity (checked above) and only be bilinear (checked in validity)
				return term.getFactor() * this->model->computeEpsilon(gamma) * getExpectationValue(term.operators[0U], gamma);
			}
			return term.getFactor() * getExpectationValue(term.operators[0U], gamma);
		};

	public:
		TermWithDOS(Utility::InputFileReader& input) : DetailModelConstructor(input) {};
	};
}