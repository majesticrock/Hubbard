#pragma once
#include "HandlerBase.hpp"

class TestHandler : public HandlerBase
{
public:
	TestHandler(mrock::Utility::InputFileReader& input, int _rank, int _numberOfRanks)
		: HandlerBase(input, _rank, _numberOfRanks) {};
	void execute(mrock::Utility::InputFileReader& input) const;
};