#pragma once
#include "DetailModelConstructor.hpp"
#include "../Models/DOSModels/BroydenDOS.hpp"
#include <complex>
#include <cmath>
#include <filesystem>
#include <mrock/utility/BinaryIO.hpp>

namespace Hubbard::Helper {
	template <class DOS>
	class TermWithDOS : protected DetailModelConstructor<Models::DOSModels::BroydenDOS<DOS>>
	{
	private:
		using Integrator = DOS::template Integrator<complex_prec>;
		Integrator _integrator{};

		// The non-summing terms (except for the epsilon terms) are proportioal to 1/#lattice sites
		// This number is probably arbitrary in the DOS-description so we will fix it to a value that works
	protected:
		// The abcissa need to be organized so that [0] -> -gamma_max + Delta gamma / 2
		// Then a spacing of Delta gamma = 2 gamma_max / Constants::BASIS_SIZE yields a proper
		// discretization that omits 0 (which is problematic when we need to compute k+Q for gamma(k) = 0
		// [n] = -gamma_max + Delta gamma * (n + 1/2)
		std::vector<global_floating_type> approximate_dos;
		const global_floating_type INV_GAMMA_DISC;

		inline global_floating_type getWeightFromIndex(size_t gamma_idx) {
			if (gamma_idx < Constants::HALF_BASIS) {
				// gamma < 0
				if (gamma_idx < Constants::EIGHTH_BASIS) {
					return this->model->getDeltaGamma();
				}
				if (gamma_idx < Constants::QUARTER_BASIS) {
					return 0.5 * this->model->getDeltaGamma();
				}
				return 0.25 * this->model->getDeltaGamma();
			}
			return getWeightFromIndex(gamma_idx - Constants::HALF_BASIS);
			//return DOS::weights_v(index);
		};

		global_floating_type getExpectationValue(const mrock::symbolic_operators::WickOperator& op, int gamma_idx) const {
			assert(op.type < mrock::symbolic_operators::OperatorType::Undefined_Type);

			int index = static_cast<int>(op.type);
			if (op.type == mrock::symbolic_operators::OperatorType::CDW_Type 
				|| op.type == mrock::symbolic_operators::OperatorType::Number_Type)
			{
				auto jt = this->wick_spin_offset.find(op.indizes[0]);
				if (jt == this->wick_spin_offset.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
				index += jt->second;
			}

			auto offset = [&op, &gamma_idx, this]() -> int {
				if (op.momentum.add_Q) {
					return this->model->shiftByQ(gamma_idx);
				}
				return gamma_idx;
				};

			//if (op.is_daggered) return std::conj(this->expecs[index](offset(), 0));
			return this->expecs[index](offset(), 0);
		};

		global_floating_type computeTerm(const mrock::symbolic_operators::WickTerm& term, int gamma_idx, int gamma_prime_idx) const {
			auto attributeVanishes = [this](int index) -> bool {
				return (!this->model->getAttributes().isFinite(index));
				};

			if (attributeVanishes(0) && attributeVanishes(1) && term.includes_type(mrock::symbolic_operators::OperatorType::CDW_Type)) {
				return global_floating_type{};
			}
			if (attributeVanishes(2) && term.includes_type(mrock::symbolic_operators::OperatorType::SC_Type)) {
				return global_floating_type{};
			}

			const global_floating_type& gamma{ this->model->getGammaFromIndex(gamma_idx) };
			const global_floating_type& gamma_prime{ this->model->getGammaFromIndex(gamma_prime_idx) };

			auto getCoefficient = [&]() {
				if (term.get_first_coefficient().depends_on('l')) {
					return gamma_prime * term.get_factor() * this->model->computeCoefficient(term.get_first_coefficient(), gamma);
				}
				return term.get_factor() * this->model->computeCoefficient(term.get_first_coefficient(), gamma);
				};

			auto getCoefficientAndExpec = [&](size_t expec_pos) {
				return getCoefficient() * getExpectationValue(term.operators[expec_pos], gamma_idx);
				};

			if (term.is_identity()) {
				if (term.has_single_coefficient()) {
					return getCoefficient();
				}
				return term.get_factor();
			}
			if (term.sums.momenta.size() > 0U) {
				if (term.is_bilinear()) {
					// bilinear
					return getCoefficient() * this->getSumOfAll(term.operators.front(), term.get_first_coefficient().depends_on('q'));
				}
				else if (term.is_quartic()) {
					// quartic
					int q_dependend = term.which_operator_depends_on('q');
					if (q_dependend < 0) throw std::invalid_argument("Suming q, but no operator depends on q.");

					return getCoefficientAndExpec(q_dependend == 0) * this->getSumOfAll(term.operators[q_dependend], term.get_first_coefficient().depends_on('q'));
				}
			}

			if (term.has_single_coefficient()) {
				// Can never be an identity (checked above) and only be bilinear or quartic (checked in validity)
				global_floating_type ret{};
				if (term.is_bilinear()) {
					const auto& op = term.operators.front();
					ret = (op.depends_on('l') ? getExpectationValue(op, gamma_prime_idx) : getExpectationValue(op, gamma_idx));
				}
				else {
					int l_dependend = term.which_operator_depends_on('l');
					ret = getExpectationValue(term.operators[l_dependend == 0], gamma_idx)
						* getExpectationValue(term.operators[l_dependend], gamma_prime_idx);
				}
				return ret * getCoefficient();
			}

			return term.get_factor() * getExpectationValue(term.operators[0U], gamma_idx);
		};

	public:
		TermWithDOS(mrock::utility::InputFileReader& input, const Models::ModelParameters& modelParameters)
			: DetailModelConstructor<Models::DOSModels::BroydenDOS<DOS>>(input, modelParameters),
			approximate_dos(Constants::BASIS_SIZE, global_floating_type{}),
			INV_GAMMA_DISC(1. / this->model->getDeltaGamma())
		{
#ifdef _EXACT_DOS
			auto dos_norm = [this]() -> global_floating_type {
				global_floating_type val{};
				for (size_t i = 0U; i < approximate_dos.size(); ++i)
				{
					val += approximate_dos[i] * getWeightFromIndex(i);
				}
				return val;
				};
#else
			auto dos_norm = [this]() -> global_floating_type {
				return this->model->getDeltaGamma() * std::reduce(approximate_dos.begin(), approximate_dos.end(), global_floating_type{});
				};
#endif

			const std::string filename = "../../data/approx_dos_dim_" + std::to_string(DOS::DIMENSION) + "_disc_" + std::to_string(Constants::BASIS_SIZE) + ".bin";
			if (std::filesystem::exists(filename)) {
				std::ifstream reader = mrock::utility::BinaryIO::readSerializedVector(approximate_dos, filename);
				if (reader.good() && std::abs(dos_norm() - 1.) < 1e-7) {
					return;
				}
				else {
					std::cerr << "An error occurred while reading the approximate dos data." << std::endl;
					std::cerr << "1 - Approximate DOS norm = " << std::scientific << std::setprecision(8) << std::abs(dos_norm() - 1.) << std::endl;
				}
			}
			else {
				std::cout << "DOS file does not exist yet. Computing..." << std::endl;
			}

			std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

			for (int i = 0; i < Constants::BASIS_SIZE; ++i)
			{
				//std::cout << i << "\t" << this->model->getGammaFromIndex(i) << "\t" << this->model->getGammaFromIndex(this->model->shiftByQ(i)) << std::endl;
				approximate_dos[i] = DOS::computeValue(this->model->getGammaFromIndex(i));
#ifdef _CLOSED_FORMULA
				if (i == 0 || i == Constants::BASIS_SIZE - 1) {
					approximate_dos[i] *= (3. / 8.);
				}
				else if (i % 3 == 0) {
					approximate_dos[i] *= (6. / 8.);
				}
				else {
					approximate_dos[i] *= (9. / 8.);
				}
#endif
			}
			global_floating_type inverse_norm = 1. / dos_norm();
			for (auto& value : approximate_dos)
			{
				value *= inverse_norm;
			}

			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << "Computed DOS for iEoM in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			if (!(mrock::utility::BinaryIO::serializeVector(approximate_dos, filename).good())) {
				std::cerr << "Error while writing dos data to " << filename << std::endl;
			}
		};
	};
}