#pragma once
#include "HandlerBase.hpp"

class PhaseHandler : public HandlerBase
{
public:
	PhaseHandler(Utility::InputFileReader& input, int _rank, int _numberOfRanks)
		: HandlerBase(input, _rank, _numberOfRanks) {};
	void execute(Utility::InputFileReader& input) const;
};