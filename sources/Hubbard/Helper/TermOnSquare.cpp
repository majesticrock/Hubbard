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
	
	complex_prec TermOnSquare::computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const
	{
		const Eigen::Vector2i l_idx = { x(l), y(l) };
		const Eigen::Vector2i k_idx = { x(k), y(k) };
		Eigen::Vector2i q_idx = { 0,0 };
		std::vector<Eigen::Vector2i> indizes = { l_idx, k_idx, q_idx };
		Eigen::Vector2i momentum_value, coeff_momentum;

		if (term.coefficients.size() > 1) throw std::invalid_argument("Undefined number of coefficients: " + std::to_string(term.coefficients.size()));
		if (term.operators.size() == 0) {
			if (term.coefficients.size() == 1) {
				coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
				return term.multiplicity * this->model->computeCoefficient(term.coefficients[0], coeff_momentum);
			}
			return static_cast<global_floating_type>(term.multiplicity);
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
				if (term.coefficients.size() == 1) {
					sumBuffer *= this->model->computeCoefficient(term.coefficients[0], coeff_momentum);
				}
				returnBuffer += sumBuffer;
			}
			return static_cast<global_floating_type>(term.multiplicity) * returnBuffer;
		};

		if (term.sum_momenta.size() > 0) {
			if (term.sum_momenta.size() > 1) throw std::invalid_argument("Too many sums: " + term.sum_momenta.size());
			if (term.operators.size() == 1) {
				// bilinear term
				if (term.sum_momenta.size() > 1) throw std::invalid_argument("Term with more than one momentum summation: " + term.sum_momenta.size());
				if (term.delta_momenta.size() == 0) throw std::invalid_argument("There is a summation without delta_kl in a bilinear term.");
				if (term.coefficients.size() == 1) {
					if (term.coefficients.back().momentum.momentum_list.size() == 0) {
						coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
						return term.multiplicity * this->model->computeCoefficient(term.coefficients[0], coeff_momentum)
							* getSumOfAll(term.operators[0]);
					}
					else {
						return compute_single_sum();
					}
				}
				return static_cast<global_floating_type>(term.multiplicity) * getSumOfAll(term.operators[0]);
			}
			if (term.operators.size() == 2) {
				// quartic term
				return compute_single_sum();
			}
			throw std::invalid_argument("There are more than 2 WickOperators: " + term.operators.size());
		}

		complex_prec returnBuffer{ 1, 0 };
		for (size_t i = 0; i < term.operators.size(); i++)
		{
			Eigen::Vector2i momentum_value = computeMomentum(term.operators[i].momentum, indizes, { 'l', 'k' });
			returnBuffer *= getExpectationValue(term.operators[i], momentum_value);
		}
		if (term.coefficients.size() == 1) {
			coeff_momentum = computeMomentum(term.coefficients[0].momentum, indizes, { 'l', 'k' });
			return term.multiplicity * this->model->computeCoefficient(term.coefficients[0], coeff_momentum) * returnBuffer;
		}
		return static_cast<global_floating_type>(term.multiplicity) * returnBuffer;
	}
}