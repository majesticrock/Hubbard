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

					for (const auto& term : wicks_M[j * number_of_basis_terms + i]) {
						checkTermValidity(term);
					}
				}
				{
					// create an input file stream and a text archive to deserialize the vector
					std::ifstream ifs(filename + "N_" + std::to_string(j + name_offset) + "_" + std::to_string(i + name_offset) + ".txt");
					boost::archive::text_iarchive ia(ifs);
					wicks_N[j * number_of_basis_terms + i].clear();
					ia >> wicks_N[j * number_of_basis_terms + i];
					ifs.close();

					for (const auto& term : wicks_N[j * number_of_basis_terms + i]) {
						checkTermValidity(term);
					}
				}
			}
		}
	}

	void ModeHelper::checkTermValidity(const SymbolicOperators::WickTerm& term)
	{
		if (term.delta_momenta.size() > 1) throw std::invalid_argument("Too many deltas: " + term.delta_momenta.size());
		if (term.delta_momenta.size() == 1) {
			if (term.delta_momenta[0].first.momentum_list.size() != 1) throw std::invalid_argument("First delta list is not of size 1: " + term.delta_momenta[0].first.momentum_list.size());
			if (term.delta_momenta[0].second.momentum_list.size() != 1) throw std::invalid_argument("To be implemented: Second delta list is not of size 1: " + term.delta_momenta[0].second.momentum_list.size());
		}

		if (term.coefficients.size() > 1U) throw std::invalid_argument("Undefined number of coefficients: " + std::to_string(term.coefficients.size()));
		if (term.operators.size() > 2U) throw std::invalid_argument("There are more than 2 WickOperators: " + term.operators.size());
		if (term.sum_momenta.size() > 0U) {
			if (!term.hasSingleCoefficient()) throw std::invalid_argument("Too many sums: " + term.sum_momenta.size());
			if (term.delta_momenta.empty()) throw std::invalid_argument("There is a summation without delta_kl.");
		}
		else {
			if (term.operators.size() > 2U) throw std::invalid_argument("A term without a sum can only be bilinear or an identity.");
		}
	}

	ModeHelper::ModeHelper(Utility::InputFileReader& input)
		: number_of_basis_terms{ input.getInt("number_of_basis_terms") }, start_basis_at{ input.getInt("start_basis_at") }
	{
		if (start_basis_at < 0) {
			// We investigate the special x-p-basis
			this->TOTAL_BASIS = Constants::BASIS_SIZE * 8;
			this->number_of_basis_terms = 10;
		}
		else {
			this->TOTAL_BASIS = Constants::BASIS_SIZE * this->number_of_basis_terms;
		}

		loadWick("../commutators/wick_");
	}
}