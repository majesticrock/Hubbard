#include "BaseModel.hpp"

namespace Hubbard
{
	BaseModel::BaseModel(const ModelParameters& _params)
		: temperature(_params.temperature), U(_params.U), V(_params.V)
	{
	}
} // namespace Hubbard