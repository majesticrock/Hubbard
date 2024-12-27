#pragma once
#include "ModeHandler.hpp"
#include "PhaseHandler.hpp"

class UnknownBoundaryHandler : public PhaseHandler, public ModeHandler
{
public:
	UnknownBoundaryHandler(mrock::Utility::InputFileReader& input, int _rank, int _numberOfRanks)
		: HandlerBase(input, _rank, _numberOfRanks), PhaseHandler(input, _rank, _numberOfRanks), ModeHandler(input, _rank, _numberOfRanks) {};
	virtual void execute(mrock::Utility::InputFileReader& input) const override;
};