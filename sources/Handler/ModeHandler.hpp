#pragma once
#include "HandlerBase.hpp"
#include "../Hubbard/Helper/ModeHelper.hpp"
#include <memory>
#include <vector>
#include <string>

class ModeHandler : virtual public HandlerBase
{
protected:
	std::unique_ptr<Hubbard::Helper::ModeHelper> getHelper(mrock::utility::InputFileReader& input, Hubbard::Models::ModelParameters& modelParameters) const;
	std::unique_ptr<Hubbard::Helper::ModeHelper> getHelper(mrock::utility::InputFileReader& input) const;

	std::vector<std::string> getFileComments(mrock::utility::InputFileReader& input, Hubbard::Helper::ModeHelper* modeHelper) const;
public:
	ModeHandler(mrock::utility::InputFileReader& input, int _rank, int _numberOfRanks)
		: HandlerBase(input, _rank, _numberOfRanks) {};
	virtual void execute(mrock::utility::InputFileReader& input) const override;
};