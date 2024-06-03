#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	void ModeHelper::checkTermValidity(const SymbolicOperators::WickTerm& term)
	{
		if (term.delta_momenta.size() > 1) throw std::invalid_argument("Too many deltas: " + term.delta_momenta.size());
		if (term.delta_momenta.size() == 1) {
			if (term.delta_momenta[0].first.momentum_list.size() != 1) throw std::invalid_argument("First delta list is not of size 1: " + term.delta_momenta[0].first.momentum_list.size());
			if (term.delta_momenta[0].second.momentum_list.size() != 1) throw std::invalid_argument("To be implemented: Second delta list is not of size 1: " + term.delta_momenta[0].second.momentum_list.size());
		}

		if (term.coefficients.size() > 1U) throw std::invalid_argument("Undefined number of coefficients: " + std::to_string(term.coefficients.size()));
		if (term.operators.size() > 2U) throw std::invalid_argument("There are more than 2 WickOperators: " + term.operators.size());
		if (term.sums.momenta.size() > 0U) {
			if (!term.hasSingleCoefficient()) throw std::invalid_argument("Too many sums: " + term.sums.momenta.size());
			if (term.delta_momenta.empty()) throw std::invalid_argument("There is a summation without delta_kl.");
		}
		else {
			if (term.operators.size() > 2U) throw std::invalid_argument("A term without a sum can only be bilinear, quartic or an identity.");
		}
		if (!(term.coefficients.empty())) {
			if (term.getFirstCoefficient().dependsOnMomentum()) {
				if (!(term.getFirstCoefficient().dependsOn('k'))) throw std::invalid_argument("Each momentum dependent term should have a k-depedance.");
			}
			if (term.getFirstCoefficient().momenta.size() > 1U)
				throw std::invalid_argument("There must not be more than 1 momentum in coefficient!");
		}
	}

	ModeHelper::ModeHelper(Utility::InputFileReader& input)
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
}