#pragma once
#include "Utility/InputFileReader.hpp"
#include "Hubbard/ModelParameters.hpp"

class HandlerBase
{
protected:
	Hubbard::ModelParameters modelParameters;
	int rank{};
	int numberOfRanks{ 1 };
public:
	HandlerBase(Utility::InputFileReader& input, int _rank, int _numberOfRanks);
	virtual ~HandlerBase() = default;
};