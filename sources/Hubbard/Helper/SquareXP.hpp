#pragma once
#include "XPModes.hpp"
#include "TermOnSquare.hpp"
#include <array>

namespace Hubbard::Helper {
	class SquareXP : public TermOnSquare, public XPModes
	{
	private:
		inline global_floating_type computeRealTerm(const SymbolicOperators::WickTerm& term, int k, int l) const {
			const auto result = this->computeTerm(term, k, l);
			if (abs(result.imag()) > ERROR_MARGIN) {
				throw std::runtime_error("computeRealTerm() encountered a complex value!");
			}
			return result.real();
		};

		virtual void fill_block_M(int i, int j) override;
		virtual void fill_block_N(int i, int j) override;

	public:
		virtual const BaseModel<global_floating_type>& getModel() const override {
			return *model;
		};
		virtual BaseModel<global_floating_type>& getModel() override {
			return *model;
		};

		SquareXP(Utility::InputFileReader& input, const ModelParameters& modelParameters) : TermOnSquare(input, modelParameters), XPModes(input)
		{};

		void setNewModelParameters(Utility::InputFileReader& input, const ModelParameters& modelParameters) override;
	};
}