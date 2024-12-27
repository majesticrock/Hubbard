#pragma once
#include "GeneralBasis.hpp"
#include "TermOnSquare.hpp"

namespace Hubbard::Helper {
	class SquareGeneral : public TermOnSquare, public GeneralBasis
	{
	private:
		virtual void fill_block_M(int i, int j) override;
		virtual void fill_block_N(int i, int j) override;

	public:
		virtual const Models::BaseModel<global_floating_type>& getModel() const override {
			return *model;
		};
		virtual Models::BaseModel<global_floating_type>& getModel() override {
			return *model;
		};

		SquareGeneral(mrock::Utility::InputFileReader& input, const Models::ModelParameters& modelParameters) : TermOnSquare(input, modelParameters), GeneralBasis(input) {};
		SquareGeneral(mrock::Utility::InputFileReader& input, std::unique_ptr<Hubbard::Models::SquareLattice::UsingBroyden>&& model_ptr)
			: TermOnSquare(std::move(model_ptr)), GeneralBasis(input) {};

		void setNewModelParameters(mrock::Utility::InputFileReader& input, const Models::ModelParameters& modelParameters) override;
	};
}