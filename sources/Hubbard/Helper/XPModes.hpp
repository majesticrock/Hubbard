#pragma once
#include "ModeHelper.hpp"

namespace Hubbard::Helper {
	class XPModes : public ModeHelper
	{
	private:
		void fillBlock(int i, int j);
	protected:
		Matrix_L K_plus, K_minus, L;
		inline double computeRealTerm(const SymbolicOperators::WickTerm& term, int l, int k) const {
			auto result = computeTerm(term, l, k);
			if (std::abs(result.imag()) > ERROR_MARGIN) {
				throw std::runtime_error("computeRealTerm() encountered a complex value!");
			}
			return result.real();
		};
		virtual void fillMatrices() override;
	public:
		XPModes(Utility::InputFileReader& input) : ModeHelper(input) { };

		virtual std::unique_ptr<std::vector<Resolvent_L>> computeCollectiveModes(std::vector<std::vector<double>>& reciever) override;
	};
}