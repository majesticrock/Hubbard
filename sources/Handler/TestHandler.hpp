#pragma once
#include "HandlerBase.hpp"

class TestHandler : public HandlerBase
{
public:
	TestHandler(mrock::utility::InputFileReader& input, int _rank, int _numberOfRanks)
		: HandlerBase(input, _rank, _numberOfRanks) {};
	void execute(mrock::utility::InputFileReader& input) const;
};