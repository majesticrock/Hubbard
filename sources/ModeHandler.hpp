#pragma once
#include "HandlerBase.hpp"

class ModeHandler : public HandlerBase
{
public:
	ModeHandler(Utility::InputFileReader& input, int _rank, int _numberOfRanks)
		: HandlerBase(input, _rank, _numberOfRanks) {};
	void execute(Utility::InputFileReader& input) const;
};