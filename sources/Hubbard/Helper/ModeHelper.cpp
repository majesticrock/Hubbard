#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	void ModeHelper::loadWick(const std::string& filename)
	{
		wicks_M.resize(number_of_basis_terms * number_of_basis_terms);
		wicks_N.resize(number_of_basis_terms * number_of_basis_terms);
		const int name_offset = (start_basis_at < 0) ? 0 : start_basis_at;

		for (int i = 0; i < number_of_basis_terms; i++)
		{
			for (int j = 0; j < number_of_basis_terms; j++)
			{
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "M_" + std::to_string(j + name_offset) + "_" + std::to_string(i + name_offset) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_M[j * number_of_basis_terms + i].clear();
					ia >> wicks_M[j * number_of_basis_terms + i];
					ifs.close();
				}
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "N_" + std::to_string(j + name_offset) + "_" + std::to_string(i + name_offset) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_N[j * number_of_basis_terms + i].clear();
					ia >> wicks_N[j * number_of_basis_terms + i];
					ifs.close();
				}
			}
		}
	}

	complex_prec ModeHelper::getExpectationValue(const SymbolicOperators::WickOperator& op, const Eigen::Vector2i& momentum_value) const
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

	complex_prec ModeHelper::getSumOfAll(const SymbolicOperators::WickOperator& op) const
	{
		auto it = wick_map.find(op.type);
		if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + op.type);

		int index = it->second;
		if (op.type == "g" || op.type == "n") {
			auto jt = wick_spin_offset.find(op.indizes[0]);
			if (jt == wick_map.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
			index += jt->second;
		}

		if (op.isDaggered) return std::conj(sum_of_all[index]);
		return sum_of_all[index];
	}

	complex_prec ModeHelper::computeTerm(const SymbolicOperators::WickTerm& term, int l, int k) const
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
			complex_prec sumBuffer = 0;
			complex_prec returnBuffer = 0;
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

		complex_prec returnBuffer = 1;
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

	ModeHelper::ModeHelper(Utility::InputFileReader& input)
	{
		std::vector<double> model_params = input.getDoubleList("model_parameters");
		Hubbard::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
			0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));

		model = std::make_unique<SquareLattice::UsingBroyden>(modelParameters);

		start_basis_at = input.getInt("start_basis_at");
		this->number_of_basis_terms = input.getInt("number_of_basis_terms");
		if (start_basis_at < 0) {
			// We investigate the special x-p-basis
			this->TOTAL_BASIS = Constants::BASIS_SIZE * 8;
			this->number_of_basis_terms = 10;
		}
		else {
			this->TOTAL_BASIS = Constants::BASIS_SIZE * this->number_of_basis_terms;
		}

		model->computePhases().print();
		model->computeExpectationValues(expecs, sum_of_all);

		loadWick("../commutators/wick_");
	}
}