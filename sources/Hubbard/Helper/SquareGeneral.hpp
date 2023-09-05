#pragma once
#include "GeneralBasis.hpp"
#include "TermOnSquare.hpp"

namespace Hubbard::Helper {
	class SquareGeneral : public GeneralBasis, TermOnSquare
	{
	private:
		virtual void fillBlock(int i, int j) override;

	public:
		virtual const BaseModel<global_floating_type>& getModel() const override {
			return *model;
		};

		SquareGeneral(Utility::InputFileReader& input) : GeneralBasis(input), TermOnSquare(input) {};
	};
}