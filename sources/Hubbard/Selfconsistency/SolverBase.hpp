#pragma once
#include "../GlobalDefinitions.hpp"
#include "../BaseModel.hpp"

namespace Hubbard::Selfconsistency {
	template <typename DataType>
	class SolverBase {
	protected:
		BaseModel<DataType>* _model{};
		ModelAttributes<DataType>* _attr{};
		const size_t NUMBER_OF_PARAMETERS;

	public:
		virtual ModelAttributes<double> computePhases(const PhaseDebuggingPolicy& debugPolicy) = 0;

		virtual ~SolverBase() = default;
		SolverBase() = delete;
		SolverBase(BaseModel<DataType>* model_ptr, ModelAttributes<DataType>* attribute_ptr)
			: _model(model_ptr), _attr(attribute_ptr), NUMBER_OF_PARAMETERS(_attr->size()) {};
	};
}