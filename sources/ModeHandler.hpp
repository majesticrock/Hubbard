#pragma once
#include "HandlerBase.hpp"
#include "Hubbard/Helper/ModeHelper.hpp"
#include <memory>
#include <vector>
#include <string>

class ModeHandler : virtual public HandlerBase
{
protected:
	std::unique_ptr<Hubbard::Helper::ModeHelper> getHelper(Utility::InputFileReader& input, Hubbard::ModelParameters& modelParameters) const;
	std::unique_ptr<Hubbard::Helper::ModeHelper> getHelper(Utility::InputFileReader& input) const;

	std::vector<std::string> getFileComments(Utility::InputFileReader& input, Hubbard::Helper::ModeHelper* modeHelper) const;
public:
	ModeHandler(Utility::InputFileReader& input, int _rank, int _numberOfRanks)
		: HandlerBase(input, _rank, _numberOfRanks) {};
	virtual void execute(Utility::InputFileReader& input) const override;
};