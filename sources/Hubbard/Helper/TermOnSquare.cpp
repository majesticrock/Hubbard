#include "TermOnSquare.hpp"

namespace Hubbard::Helper {
	complex_prec TermOnSquare::getExpectationValue(const SymbolicOperators::WickOperator& op, const Eigen::Vector2i& momentum_value) const
	{
		auto it = wick_map.find(op.type);
		if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + op.type);

		int index = it->second;
		if (op.type == "g" || op.type == "n") {
			auto jt = wick_spin_offset.find(op.indizes[0]);
			if (jt == wick_map.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
			index += jt->second;
		}

		if (op.isDaggered) return std::conj(expecs[index](momentum_value(0), momentum_value(1)));
		return expecs[index](momentum_value(0), momentum_value(1));
	}

	Eigen::Vector2i TermOnSquare::computeMomentum(const SymbolicOperators::Momentum& momentum, const std::vector<Eigen::Vector2i>& indizes, const std::vector<char>& momenta) const {
		Eigen::Vector2i buffer = { 0,0 };
		for (int i = 0; i < momenta.size(); ++i)
		{
			int mom_idx = momentum.isUsed(momenta[i]);
			if (mom_idx < 0) continue;
			buffer += momentum.momentum_list[mom_idx].first * indizes[i];
		}
		if (momentum.add_Q) {
			buffer(0) += Constants::K_DISCRETIZATION;
			buffer(1) += Constants::K_DISCRETIZATION;
		}
		clean_factor_2pi(buffer);
		return buffer;
	}

	complex_prec TermOnSquare::computeTerm(const SymbolicOperators::WickTerm& term, int k, int l) const
	{
		const Eigen::Vector2i l_idx = { x(l), y(l) };
		const Eigen::Vector2i k_idx = { x(k), y(k) };
		Eigen::Vector2i q_idx = { 0,0 };
		std::vector<Eigen::Vector2i> indizes = { l_idx, k_idx, q_idx };
		Eigen::Vector2i momentum_value, coeff_momentum;

		if (term.isIdentity()) {
			if (term.hasSingleCoefficient()) {
				coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
				return term.getFactor() * this->model->computeCoefficient(term.coefficients[0], coeff_momentum);
			}
			return term.getFactor();
		}

		auto compute_single_sum = [&]() -> complex_prec {
			complex_prec sumBuffer{};
			complex_prec returnBuffer{};
			for (int q = 0; q < Constants::BASIS_SIZE; q++)
			{
				q_idx = { x(q), y(q) };
				indizes = { l_idx, k_idx, q_idx };
				sumBuffer = 1;
				for (size_t i = 0; i < term.operators.size(); i++)
				{
					momentum_value = computeMomentum(term.operators[i].momentum, indizes, { 'l', 'k', 'q' });
					sumBuffer *= getExpectationValue(term.operators[i], momentum_value);
				}
				coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k', 'q' });
				if (term.hasSingleCoefficient()) {
					sumBuffer *= this->model->computeCoefficient(term.coefficients[0], coeff_momentum);
				}
				returnBuffer += sumBuffer;
			}
			return term.getFactor() * returnBuffer;
			};

		if (term.sum_momenta.size() > 0U) {
			if (term.isBilinear()) {
				// bilinear term
				if (term.hasSingleCoefficient()) {
					if (term.coefficients.back().dependsOnMomentum()) {
						return compute_single_sum();
					}
					else {
						coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
						return term.getFactor() * this->model->computeCoefficient(term.coefficients[0], coeff_momentum) * getSumOfAll(term.operators[0]);
					}
				}
				return term.getFactor() * getSumOfAll(term.operators[0]);
			}
			if (term.isQuartic()) {
				// quartic term
				return compute_single_sum();
			}
		}

		complex_prec returnBuffer{ 1, 0 };
		for (size_t i = 0U; i < term.operators.size(); ++i)
		{
			Eigen::Vector2i momentum_value = computeMomentum(term.operators[i].momentum, indizes, { 'l', 'k' });
			returnBuffer *= getExpectationValue(term.operators[i], momentum_value);
		}
		if (term.hasSingleCoefficient()) {
			coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
			return term.getFactor() * this->model->computeCoefficient(term.coefficients[0], coeff_momentum) * returnBuffer;
		}
		return term.getFactor() * returnBuffer;
	}
}