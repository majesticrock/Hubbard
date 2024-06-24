#pragma once
#include <Utility/InputFileReader.hpp>
#include "../Hubbard/Models/ModelParameters.hpp"
#include <string>

class HandlerBase
{
protected:
	Hubbard::Models::ModelParameters modelParameters;
	int rank{};
	int numberOfRanks{ 1 };

	inline std::string getOutputFolder(Utility::InputFileReader& input) const {
		return input.getString("lattice_type") + "/" + input.getString("output_folder");
	}
public:
	HandlerBase(Utility::InputFileReader& input, int _rank, int _numberOfRanks);
	virtual ~HandlerBase() = default;
	virtual void execute(Utility::InputFileReader& input) const = 0;
};