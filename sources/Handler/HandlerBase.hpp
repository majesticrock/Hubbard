#pragma once
#include <mrock/utility/InputFileReader.hpp>
#include "../Hubbard/Models/ModelParameters.hpp"
#include <string>

class HandlerBase
{
protected:
	Hubbard::Models::ModelParameters modelParameters;
	int rank{};
	int numberOfRanks{ 1 };

	inline std::string getOutputFolder(mrock::utility::InputFileReader& input) const {
		return input.getString("lattice_type") + "/" + input.getString("output_folder");
	}
public:
	HandlerBase(mrock::utility::InputFileReader& input, int _rank, int _numberOfRanks);
	virtual ~HandlerBase() = default;
	virtual void execute(mrock::utility::InputFileReader& input) const = 0;
};