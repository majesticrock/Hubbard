#include "TermOnSquare.hpp"

namespace Hubbard::Helper {
	global_floating_type TermOnSquare::getExpectationValue(const mrock::symbolic_operators::WickOperator& op, const Eigen::Vector2i& momentum_value) const
	{
		assert(op.type < mrock::symbolic_operators::OperatorType::Undefined_Type);

		int index = static_cast<int>(op.type);
		if (op.type == mrock::symbolic_operators::OperatorType::CDW_Type 
			|| op.type == mrock::symbolic_operators::OperatorType::Number_Type) 
		{
			auto jt = wick_spin_offset.find(op.indizes[0]);
			if (jt == wick_spin_offset.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
			index += jt->second;
		}

		if (op.is_daggered) return expecs[index](momentum_value(0), momentum_value(1));
		return expecs[index](momentum_value(0), momentum_value(1));
	}

	Eigen::Vector2i TermOnSquare::compute_momentum_list(const mrock::symbolic_operators::MomentumList& momentum, const Eigen::Vector2i& k, const Eigen::Vector2i& l) const
	{
		if (momentum.empty()) {
			return Eigen::Vector2i::Zero();
		}
		return compute_momentum_no_q(momentum.front(), k, l);
	}

	Eigen::Vector2i TermOnSquare::compute_momentum_no_q(const mrock::symbolic_operators::Momentum& momentum, const Eigen::Vector2i& k, const Eigen::Vector2i& l) const
	{
		Eigen::Vector2i momentum_value = Eigen::Vector2i::Zero();
		for (const auto& momentum_symbol : momentum.momentum_list) {
			switch(static_cast<char>(momentum_symbol.name)) {
			case 'k':
				momentum_value += momentum_symbol.factor * k;
				break;
			case 'l':
				momentum_value += momentum_symbol.factor * l;
				break;
			case 'x':
				momentum_value += momentum_symbol.factor * this->mode_momentum;
				break;
			default:
				break;
			}
		}
		if(momentum.add_Q) {
			momentum_value(0) += Constants::K_DISCRETIZATION;
			momentum_value(1) += Constants::K_DISCRETIZATION;
		}
		clean_factor_2pi(momentum_value);
		return momentum_value;
	}

	global_floating_type TermOnSquare::compute_sum(const mrock::symbolic_operators::WickTerm& term, const Eigen::Vector2i& k, const Eigen::Vector2i& l) const {
		global_floating_type value{ static_cast<global_floating_type>(term.multiplicity) };

		// Compute coefficient
		assert(term.has_single_coefficient());
		int sum_of_all_index{ 0 };
		mrock::symbolic_operators::Coefficient const& coeff = term.coefficients.front();
		if(coeff.momenta.empty()) {
			value *= this->model->compute_coefficient(coeff, Eigen::Vector2i::Zero());
		}
		else {
			value *= this->model->compute_coefficient(coeff, compute_momentum_no_q(coeff.momenta.front(), k, l));
			if(coeff.depends_on('q')) sum_of_all_index = 1;
		}

		// Compute sum
		const int q_dependend = term.which_operator_depends_on('q');
		mrock::symbolic_operators::WickOperator const* const summed_op = &(term.operators[q_dependend]);
		mrock::symbolic_operators::WickOperator const* const other_op = term.is_bilinear() ? nullptr : &(term.operators[q_dependend == 0]);
		int index = static_cast<int>(summed_op->type);
		if (summed_op->type == mrock::symbolic_operators::OperatorType::CDW_Type 
			|| summed_op->type == mrock::symbolic_operators::OperatorType::Number_Type) 
		{
			auto jt = wick_spin_offset.find(summed_op->indizes[0]);
			if (jt == wick_spin_offset.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
			index += jt->second;
		}
		value *= sum_of_all(index, sum_of_all_index);

		// Compute the expectation value that is not summed over, if it exists
		if (other_op) {
			value *= this->getExpectationValue(*other_op, compute_momentum_no_q(other_op->momentum, k, l));
		}
		return value;
	}

	global_floating_type TermOnSquare::computeTerm(const mrock::symbolic_operators::WickTerm& term, const int k, const int l) const
	{
		const OrderType model_order = this->model->get_order();
		if(term.includes_type(mrock::symbolic_operators::OperatorType::CDW_Type) && !(model_order & (OrderType::CDW | OrderType::AFM))) return global_floating_type{};
		if(term.includes_type(mrock::symbolic_operators::OperatorType::SC_Type)  && !(model_order & OrderType::SC)) return global_floating_type{};

		const std::vector<char>& momenta_plain = { 'k', 'l' };
		std::vector<Eigen::Vector2i> indizes = { { x(k), y(k) }, { x(l), y(l) } };
		Eigen::Vector2i coeff_momentum;

		if (term.is_identity()) {
			if (term.has_single_coefficient()) {
				coeff_momentum = compute_momentum_list(term.coefficients[0].momenta, indizes[0], indizes[1]);
				return term.get_factor() * this->model->compute_coefficient(term.coefficients[0], coeff_momentum);
			}
			return term.get_factor();
		}

		if (term.sums.has_momentum()) {
			return compute_sum(term, indizes[0], indizes[1]);
		}

		global_floating_type returnBuffer{ 1 };
		for (const auto& op : term.operators)
		{
			returnBuffer *= getExpectationValue(op, compute_momentum_no_q(op.momentum, indizes[0], indizes[1]));
		}
		if (term.has_single_coefficient()) {
			coeff_momentum = compute_momentum_list(term.coefficients[0].momenta, indizes[0], indizes[1]);
			return term.get_factor() * this->model->compute_coefficient(term.coefficients[0], coeff_momentum) * returnBuffer;
		}
		return term.get_factor() * returnBuffer;
	}
}