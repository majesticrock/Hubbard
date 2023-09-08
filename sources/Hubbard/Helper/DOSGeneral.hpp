#pragma once
#include "GeneralBasis.hpp"
#include "TermWithDOS.hpp"

namespace Hubbard::Helper {
	class DOSGeneral : public GeneralBasis, public TermWithDOS
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

		DOSGeneral(Utility::InputFileReader& input) : GeneralBasis(input), TermWithDOS(input) {};
	};
}