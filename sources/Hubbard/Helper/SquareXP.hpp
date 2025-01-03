#pragma once
#include "XPModes.hpp"
#include "TermOnSquare.hpp"
#include <array>

namespace Hubbard::Helper {
	class SquareXP : public TermOnSquare, public XPModes
	{
	private:
		inline global_floating_type computeRealTerm(const mrock::symbolic_operators::WickTerm& term, int k, int l) const {
			const auto result = this->computeTerm(term, k, l);
			if (abs(std::imag(result)) > ERROR_MARGIN) {
				throw std::runtime_error("computeRealTerm() encountered a complex value!");
			}
			return std::real(result);
		};

		virtual void fill_block_M(int i, int j) override;
		virtual void fill_block_N(int i, int j) override;

	public:
		virtual const Models::BaseModel<global_floating_type>& getModel() const override {
			return *model;
		};
		virtual Models::BaseModel<global_floating_type>& getModel() override {
			return *model;
		};

		SquareXP(mrock::utility::InputFileReader& input, const Models::ModelParameters& modelParameters, const Eigen::Vector2i& _mode_momentum = { 0, 0 })
			: TermOnSquare(input, modelParameters, _mode_momentum), XPModes(input)
		{};

		void setNewModelParameters(mrock::utility::InputFileReader& input, const Models::ModelParameters& modelParameters) override;
	};
}