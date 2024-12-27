#pragma once
#include <mrock/Utility/InputFileReader.hpp>
#include "../Hubbard/Models/ModelParameters.hpp"
#include <string>

class HandlerBase
{
protected:
	Hubbard::Models::ModelParameters modelParameters;
	int rank{};
	int numberOfRanks{ 1 };

	inline std::string getOutputFolder(mrock::Utility::InputFileReader& input) const {
		return input.getString("lattice_type") + "/" + input.getString("output_folder");
	}
public:
	HandlerBase(mrock::Utility::InputFileReader& input, int _rank, int _numberOfRanks);
	virtual ~HandlerBase() = default;
	virtual void execute(mrock::Utility::InputFileReader& input) const = 0;
};