#pragma once
#include "GeneralBasis.hpp"
#include "TermOnSquare.hpp"

namespace Hubbard::Helper {
	class SquareGeneral : public GeneralBasis, public TermOnSquare
	{
	private:
		virtual void fillBlock(int i, int j) override;

	public:
		virtual const BaseModel<global_floating_type>& getModel() const override {
			return *model;
		};
		virtual BaseModel<global_floating_type>& getModel() override {
			return *model;
		};

		SquareGeneral(Utility::InputFileReader& input) : GeneralBasis(input), TermOnSquare(input) {};
		SquareGeneral(Utility::InputFileReader& input, std::unique_ptr<Hubbard::SquareLattice::UsingBroyden>&& model_ptr) 
			: GeneralBasis(input), TermOnSquare(std::move(model_ptr)) {};
	};
}