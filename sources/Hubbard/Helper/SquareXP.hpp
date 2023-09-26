#pragma once
#include "XPModes.hpp"
#include "TermOnSquare.hpp"

namespace Hubbard::Helper {
	class SquareXP : public XPModes, public TermOnSquare
	{
	private:
		inline global_floating_type computeRealTerm(const SymbolicOperators::WickTerm& term, int l, int k) const {
			auto result = computeTerm(term, l, k);
			if (abs(result.imag()) > ERROR_MARGIN) {
				throw std::runtime_error("computeRealTerm() encountered a complex value!");
			}
			return result.real();
		};

		virtual void fillBlock(int i, int j) override;

	public:
		virtual const BaseModel<global_floating_type>& getModel() const override {
			return *model;
		};
		virtual BaseModel<global_floating_type>& getModel() override {
			return *model;
		};

		SquareXP(Utility::InputFileReader& input) : XPModes(input), TermOnSquare(input) {};
	};
}