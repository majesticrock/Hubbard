#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	void ModeHelper::checkTermValidity(const mrock::symbolic_operators::WickTerm& term)
	{
		if (term.delta_momenta.size() > 1) throw mrock::symbolic_operators::bad_term_exception("Too many deltas: " + term.delta_momenta.size(), term);
		if (term.delta_momenta.size() == 1) {
			if (term.delta_momenta[0].first.momentum_list.size() != 1) throw mrock::symbolic_operators::bad_term_exception("First delta list is not of size 1: " + term.delta_momenta[0].first.momentum_list.size(), term);
		}

		if (term.coefficients.size() > 1U) throw mrock::symbolic_operators::bad_term_exception("Undefined number of coefficients: " + std::to_string(term.coefficients.size()), term);
		if (term.operators.size() > 2U) throw mrock::symbolic_operators::bad_term_exception("There are more than 2 WickOperators: " + term.operators.size(), term);
		if (term.sums.momenta.size() > 0U) {
			if (!term.has_single_coefficient()) throw mrock::symbolic_operators::bad_term_exception("Too many sums: " + term.sums.momenta.size(), term);
			if (term.delta_momenta.empty()) throw mrock::symbolic_operators::bad_term_exception("There is a summation without delta_kl.", term);
		}
		else {
			if (term.operators.size() > 2U) throw mrock::symbolic_operators::bad_term_exception("A term without a sum can only be bilinear, quartic or an identity.", term);
		}
		if (!(term.coefficients.empty())) {
			//if (term.get_first_coefficient().depends_on_momentum()) {
			//	if (!(term.get_first_coefficient().depends_on('k'))) throw mrock::symbolic_operators::bad_term_exception("Each momentum dependent term should have a k-depedance.", term);
			//}
			if (term.get_first_coefficient().momenta.size() > 1U)
				throw mrock::symbolic_operators::bad_term_exception("There must not be more than 1 momentum in coefficient!", term);
		}
	}

	ModeHelper::ModeHelper(mrock::utility::InputFileReader& input)
		: number_of_basis_terms{ input.getInt("number_of_basis_terms") }, start_basis_at{ input.getInt("start_basis_at") }, usingDOS(input.getBool("use_DOS"))
	{
		if (this->usingDOS) {
			const auto lattice = input.getString("lattice_type");
			if (lattice == "square") {
				this->dos_dimension = 2;
			}
			else if (lattice == "cube") {
				this->dos_dimension = 3;
			}
			else {
				std::cerr << "Did not recognize lattice type in ModeHelper!" << std::endl;
				throw;
			}
		}
		if (this->start_basis_at < 0) {
			// We investigate the special x-p-basis
			this->TOTAL_BASIS = Constants::BASIS_SIZE * 10;
			this->number_of_basis_terms = 12;
		}
		else {
			this->TOTAL_BASIS = Constants::BASIS_SIZE * this->number_of_basis_terms;
		}

		const std::string prefix = this->start_basis_at < 0 ? "XP_" : "";
		const std::string subfolder = input.getString("compute_what") == "dispersions" ? "dispersions/" : "";
		wicks.load("../commutators/hubbard/" + subfolder, this->start_basis_at < 0, this->number_of_basis_terms, this->start_basis_at);
		try {
			for (const auto& collector : wicks.M) {
				for (const auto& term : collector) {
					this->checkTermValidity(term);
				}
			}
			for (const auto& collector : wicks.N) {
				for (const auto& term : collector) {
					this->checkTermValidity(term);
				}
			}
		}
		catch (mrock::symbolic_operators::bad_term_exception const& ex) {
			std::cout << "Encountered exception on term " << ex.which_term() << std::endl;
			throw ex;
		}
	}
}