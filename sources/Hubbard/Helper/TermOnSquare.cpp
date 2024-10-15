#include "TermOnSquare.hpp"

namespace Hubbard::Helper {
	global_floating_type TermOnSquare::getExpectationValue(const SymbolicOperators::WickOperator& op, const Eigen::Vector2i& momentum_value) const
	{
		assert(op.type < SymbolicOperators::Undefined_Type);

		int index = static_cast<int>(op.type);
		if (op.type == SymbolicOperators::CDW_Type || op.type == SymbolicOperators::Number_Type) {
			auto jt = wick_spin_offset.find(op.indizes[0]);
			if (jt == wick_spin_offset.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
			index += jt->second;
		}

		if (op.isDaggered) return expecs[index](momentum_value(0), momentum_value(1));
		return expecs[index](momentum_value(0), momentum_value(1));
	}

	Eigen::Vector2i TermOnSquare::computeMomentum(const SymbolicOperators::MomentumList& momentum,
		const std::vector<Eigen::Vector2i>& indizes, const std::vector<char>& momenta) const
	{
		if (momentum.empty()) {
			return Eigen::Vector2i::Zero();
		}
		const SymbolicOperators::Momentum& mom = momentum.front();
		int mom_idx = mom.isUsed('x');
		Eigen::Vector2i buffer = mom_idx > 0 ? (mom.momentum_list[mom_idx].first * this->mode_momentum).eval() : Eigen::Vector2i::Zero();
		for (int i = 0; i < momenta.size(); ++i)
		{
			int mom_idx = mom.isUsed(momenta[i]);
			if (mom_idx < 0) continue;
			buffer += mom.momentum_list[mom_idx].first * indizes[i];
		}
		if (mom.add_Q) {
			buffer(0) += Constants::K_DISCRETIZATION;
			buffer(1) += Constants::K_DISCRETIZATION;
		}
		clean_factor_2pi(buffer);
		return buffer;
	}

	global_floating_type TermOnSquare::compute_simple_sum(const SymbolicOperators::WickTerm& term, int k, int l) const {
		const int q_dependend = term.whichOperatorDependsOn('q');
		SymbolicOperators::WickOperator const* const summed_op = &(term.operators[q_dependend]);
		SymbolicOperators::WickOperator const* const other_op = term.isBilinear() ? nullptr : &(term.operators[q_dependend == 0]);

		int index = static_cast<int>(summed_op->type);
		if (summed_op->type == SymbolicOperators::CDW_Type || summed_op->type == SymbolicOperators::Number_Type) {
			auto jt = wick_spin_offset.find(op.indizes[0]);
			if (jt == wick_spin_offset.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
			index += jt->second;
		}
		global_floating_type value{ term.multiplicity * sum_of_all(index, 0) };
		
	}

	global_floating_type TermOnSquare::computeTerm(const SymbolicOperators::WickTerm& term, int k, int l) const
	{
		const std::vector<char>& momenta_plain = { 'k', 'l' };
		const std::vector<char>& momenta_summed = { 'k', 'l', 'q' };
		std::vector<Eigen::Vector2i> indizes = { { x(k), y(k) }, { x(l), y(l) }, { 0, 0 } };
		Eigen::Vector2i momentum_value, coeff_momentum;

		if (term.isIdentity()) {
			if (term.hasSingleCoefficient()) {
				coeff_momentum = computeMomentum(term.coefficients[0].momenta, indizes, momenta_plain);
				return term.getFactor() * this->model->computeCoefficient(term.coefficients[0], coeff_momentum);
			}
			return term.getFactor();
		}

		auto compute_single_sum = [&]() -> global_floating_type {
			global_floating_type sumBuffer{};
			global_floating_type returnBuffer{};
			for (int q = 0; q < Constants::BASIS_SIZE; q++)
			{
				indizes[2] = { x(q), y(q) };
				sumBuffer = 1;
				for (size_t i = 0; i < term.operators.size(); i++)
				{
					momentum_value = computeMomentum(SymbolicOperators::MomentumList(term.operators[i].momentum), indizes, momenta_summed);
					sumBuffer *= getExpectationValue(term.operators[i], momentum_value);
				}
				coeff_momentum = computeMomentum(term.coefficients[0].momenta, indizes, momenta_summed);
				if (term.hasSingleCoefficient()) {
					sumBuffer *= this->model->computeCoefficient(term.coefficients[0], coeff_momentum);
				}
				returnBuffer += sumBuffer;
			}
			return term.getFactor() * returnBuffer;
			};

		if (term.sums.momenta.size() > 0U) {
			if (term.isBilinear()) {
				// bilinear term
				if (term.hasSingleCoefficient()) {
					if (term.coefficients.back().dependsOnMomentum()) {
						return compute_single_sum();
					}
					else {
						coeff_momentum = computeMomentum(term.coefficients[0].momenta, indizes, momenta_plain);
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

		global_floating_type returnBuffer{ 1 };
		for (size_t i = 0U; i < term.operators.size(); ++i)
		{
			Eigen::Vector2i momentum_value = computeMomentum(SymbolicOperators::MomentumList(term.operators[i].momentum), indizes, momenta_plain);
			returnBuffer *= getExpectationValue(term.operators[i], momentum_value);
		}
		if (term.hasSingleCoefficient()) {
			coeff_momentum = computeMomentum(term.coefficients[0].momenta, indizes, momenta_plain);
			return term.getFactor() * this->model->computeCoefficient(term.coefficients[0], coeff_momentum) * returnBuffer;
		}
		return term.getFactor() * returnBuffer;
	}
}