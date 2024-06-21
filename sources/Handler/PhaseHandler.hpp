#pragma once
#include "HandlerBase.hpp"
#include "../Hubbard/Constants.hpp"

class PhaseHandler : virtual public HandlerBase
{
protected:
	std::vector<double> model_params;
	double GLOBAL_IT_LIMS[2] = { 0, 0 };
	double FIRST_IT_RANGE{};
	double FIRST_IT_MIN{};
	double FIRST_IT_MAX{};
	const int GLOBAL_IT_STEPS{};
	const int FIRST_IT_STEPS{};

public:
	PhaseHandler(Utility::InputFileReader& input, int _rank, int _numberOfRanks)
		: HandlerBase(input, _rank, _numberOfRanks), model_params{ input.getDoubleList("model_parameters") },
		GLOBAL_IT_STEPS{ input.getInt("global_iterator_steps") }, FIRST_IT_STEPS{ GLOBAL_IT_STEPS / _numberOfRanks }
	{
		// Setup the number of steps
		GLOBAL_IT_LIMS[1] = input.getDouble("global_iterator_upper_limit");
		for (int i = 0; i < Hubbard::Constants::option_list.size(); ++i)
		{
			if (input.getString("global_iterator_type") == Hubbard::Constants::option_list[i]) {
				GLOBAL_IT_LIMS[0] = model_params[i];
				FIRST_IT_RANGE = (GLOBAL_IT_LIMS[1] - GLOBAL_IT_LIMS[0]) / numberOfRanks;
				FIRST_IT_MIN = GLOBAL_IT_LIMS[0] + rank * FIRST_IT_RANGE;
				FIRST_IT_MAX = FIRST_IT_MIN + FIRST_IT_RANGE;
				model_params[i] = FIRST_IT_MIN;
			}
		}
	};
	virtual void execute(Utility::InputFileReader& input) const override;
};