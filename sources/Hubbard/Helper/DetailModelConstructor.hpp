#pragma once
#include <mrock/utility/InputFileReader.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>
#include <map>
#include <memory>
#include "../Models/ModelParameters.hpp"
#include "../Models/ModelAttributes.hpp"

namespace Hubbard::Helper {
	namespace DetailModelConstructorSettings {
		inline bool print_mean_field_result = true;
	}

	template <class Model>
	class DetailModelConstructor {
	private:
		void initialize_model(mrock::utility::InputFileReader& input, const Models::ModelParameters& modelParameters) {
			if (input.getString("ratio_CDW_SC") != "-1") {
				model->set_CDW_SC_ratio(input.getDouble("ratio_CDW_SC"));
			}
			Models::ModelAttributes<global_floating_type> result = model->computePhases();

			if (modelParameters.U > 0 && modelParameters.V > 0) { // only for U>0 and V>0
				if (result.isFinite(0) || result.isFinite(1)) {
					decltype(result) copy{ result };
					if (result.isFinite(0)) {
						copy[1] = result[0];
						copy[0] = 0;
					}
					else {
						copy[0] = result[1];
						copy[1] = 0;
					}

					decltype(model) model_copy_ptr = std::make_unique<Model>(modelParameters, copy);
					copy = model_copy_ptr->computePhases();
					if (copy.converged) {
						if (model_copy_ptr->freeEnergyPerSite() <= model->freeEnergyPerSite()) {
							model = std::move(model_copy_ptr);
						}
					}
				}
			}
			if (DetailModelConstructorSettings::print_mean_field_result)
			{
				std::cout << "Solution for " << model->parametersAsTriplet() << std::endl;
				model->getAttributes().print();
			}
			model->computeExpectationValues(expecs, sum_of_all);
		};
	protected:
		std::vector<ValueArray> expecs{};
		ValueArray sum_of_all;
		std::unique_ptr<Model> model{};

		const std::map<mrock::symbolic_operators::Index, int> wick_spin_offset = {
			{mrock::symbolic_operators::Index::SpinUp, 0}, 
			{mrock::symbolic_operators::Index::SpinDown, 4}, 
			{mrock::symbolic_operators::Index::Sigma, 6} };

		inline global_floating_type getSumOfAll(const mrock::symbolic_operators::WickOperator& op, int cos_modulation = 0) const {
			assert(op.type < mrock::symbolic_operators::OperatorType::Undefined_Type);

			int index = static_cast<int>(op.type);
			if (op.type == mrock::symbolic_operators::OperatorType::CDW_Type 
				|| op.type == mrock::symbolic_operators::OperatorType::Number_Type) 
			{
				auto jt = wick_spin_offset.find(op.indizes[0]);
				if (jt == wick_spin_offset.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
				index += jt->second;
			}

			if (op.is_daggered) return sum_of_all(index, cos_modulation);
			return sum_of_all(index, cos_modulation);
		};
		void internal_setNewModelParameters(mrock::utility::InputFileReader& input, const Models::ModelParameters& modelParameters)
		{
			this->model->setNewModelParameters(modelParameters, Model::SYSTEM_TYPE);
			this->initialize_model(input, modelParameters);
		};

	public:
		virtual ~DetailModelConstructor() = default;

		DetailModelConstructor(mrock::utility::InputFileReader& input, const Models::ModelParameters& modelParameters)
			: model(std::make_unique<Model>(modelParameters))
		{
			this->initialize_model(input, modelParameters);
		};

		DetailModelConstructor(std::unique_ptr<Model>&& model_ptr) : model{ std::move(model_ptr) }
		{
			model->computeExpectationValues(expecs, sum_of_all);
		};
	};
}