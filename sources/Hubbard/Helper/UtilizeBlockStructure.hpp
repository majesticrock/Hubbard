#pragma once

#include "DOS_XP.hpp"
#include "../DensityOfStates/Square.hpp"

namespace Hubbard::Helper {
	using ParentClass = DOS_XP<DensityOfStates::Square>;
	class UtilizeBlockStructure : public ParentClass {

	public:
		UtilizeBlockStructure(Utility::InputFileReader& input, const ModelParameters& modelParameters)
			: ParentClass(input, modelParameters) 
		{};

		virtual std::vector<ResolventReturnData> computeCollectiveModes(std::vector<std::vector<global_floating_type>>& reciever) override;
	};
}