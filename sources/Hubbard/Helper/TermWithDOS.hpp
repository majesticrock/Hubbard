#pragma once
#include "DetailModelConstructor.hpp"
#include "../DOSModels/BroydenDOS.hpp"
#include <complex>

namespace Hubbard::Helper {
	template <class DOS>
	class TermWithDOS : protected DetailModelConstructor<DOSModels::BroydenDOS<DOS>>
	{
	private:
		using Integrator = typename DOS::Integrator<complex_prec>;
		Integrator _integrator{};

	protected:
		complex_prec getExpectationValue(const SymbolicOperators::WickOperator& op, int gamma_idx) const {
			auto it = this->wick_map.find(op.type);
			if (it == this->wick_map.end()) throw std::invalid_argument("Term type not recognized: " + op.type);

			int index = it->second;
			if (op.type == "g" || op.type == "n") {
				auto jt = this->wick_spin_offset.find(op.indizes[0]);
				if (jt == this->wick_map.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
				index += jt->second;
			}

			auto offset = [&op, &gamma_idx]() -> int {
				if(op.momentum.add_Q){
					if(gamma_idx < DOS::size()){
						return DOS::size();
					}
					return -DOS::size();
				}
				return 0;
			};

			if (op.isDaggered) return std::conj(this->expecs[index](gamma_idx + offset(), 0));
			return this->expecs[index](gamma_idx + offset(), 0);
		};

		complex_prec computeTerm(const SymbolicOperators::WickTerm& term, int gamma_idx, int gamma_prime_idx) {
			const global_floating_type gamma{ DOS::abscissa_v(gamma_idx) };
			// For now terms that include l (by extend gamma') are excluded.
			// We expect these terms to be 0.
			//const global_floating_type gamma_prime{ DOS::abscissa_v(gamma_prime_idx) };

			if (!(term.coefficients.empty()) && term.sum_momenta.empty()) {
				if (term.coefficients.front().name != "\\epsilon_0") {
					// This we should be able to choose this as zero?
					// To be verified; idea is that these terms scale with 1/N and N can be set
					// arbitrarily in the DOS-based formalism
					return 0;
				}
			}

			auto getCoefficient = [&](global_floating_type x) {
				if(term.coefficients.back().momentum.add_Q){
					return this->model->computeCoefficient(term.coefficients.back(), -x);
				}
				return this->model->computeCoefficient(term.coefficients.back(), x);
			};

			auto getCoefficientAndExpec = [&](int x_idx, size_t expec_pos) {
				return this->model->computeCoefficient(term.coefficients.back(), DOS::abscissa_v(x_idx)) 
					* getExpectationValue(term.operators[expec_pos], x_idx);
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
								* (getExpectationValue(term.operators[0U], index)
									- getExpectationValue(term.operators[0U], index + DOS::size()));
						};
						return term.getFactor() * getCoefficient(gamma) * _integrator.integrate_by_index(_integrate_lambda);
					}
					return term.getFactor() * getCoefficient(gamma) * this->getSumOfAll(term.operators.front());
				}
				else if (term.isQuartic()) {
					// quartic
					int q_dependend = term.whichOperatorDependsOn('q');
					if (q_dependend < 0) throw std::invalid_argument("Suming q, but no operator depends on q.");

					if (term.coefficients.back().dependsOn('q')) {
						auto _integrate_lambda = [&](size_t index) -> complex_prec {
							return DOS::abscissa_v(index)
								* (getExpectationValue(term.operators[q_dependend], index)
									- getExpectationValue(term.operators[q_dependend], index + DOS::size()));
						};
						// q_dependend can either be 1 or 0; the other operator always depends solely on k
						// Hence q_dependend == 0 gives the positions of the k-dependend operator
						return term.getFactor() * getCoefficientAndExpec(gamma_idx, q_dependend == 0)
							* _integrator.integrate_by_index(_integrate_lambda);
					}
					return term.getFactor()	* getCoefficientAndExpec(gamma_idx, q_dependend == 0) * this->getSumOfAll(term.operators[q_dependend]);
				}
			}

			if (term.hasSingleCoefficient()) {
				if (term.coefficients.back().name != "\\epsilon_0") {
					// This we should be able to choose this as zero?
					// To be verified
					return 0;
				}
				// Can never be an identity (checked above) and only be bilinear (checked in validity)
				return term.getFactor() * getCoefficient(gamma) * getExpectationValue(term.operators[0U], gamma_idx);
			}
			return term.getFactor() * getExpectationValue(term.operators[0U], gamma_idx);
		};

	public:
		TermWithDOS(Utility::InputFileReader& input) : DetailModelConstructor<DOSModels::BroydenDOS<DOS>>(input) {};
	};
}