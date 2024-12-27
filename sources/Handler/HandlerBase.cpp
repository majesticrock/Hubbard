#include "HandlerBase.hpp"
#include <vector>

HandlerBase::HandlerBase(mrock::Utility::InputFileReader& input, int _rank, int _numberOfRanks)
	: modelParameters(input.getDoubleList("model_parameters"), 0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type")),
	rank(_rank), numberOfRanks(_numberOfRanks)
{
}