#pragma once
#include <Utility/InputFileReader.hpp>
#include <SymbolicOperators/WickTerm.hpp>
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
		void initialize_model(Utility::InputFileReader& input, const Models::ModelParameters& modelParameters) {
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

		const std::map<SymbolicOperators::Index, int> wick_spin_offset = {
			{SymbolicOperators::SpinUp, 0}, {SymbolicOperators::SpinDown, 4}, {SymbolicOperators::Sigma, 6} };

		inline complex_prec getSumOfAll(const SymbolicOperators::WickOperator& op, int cos_modulation = 0) const {
			assert(op.type < SymbolicOperators::Undefined_Type);

			int index = static_cast<int>(op.type);
			if (op.type == SymbolicOperators::CDW_Type || op.type == SymbolicOperators::Number_Type) {
				auto jt = wick_spin_offset.find(op.indizes[0]);
				if (jt == wick_spin_offset.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
				index += jt->second;
			}

			if (op.isDaggered) return std::conj(sum_of_all(index, cos_modulation));
			return sum_of_all(index, cos_modulation);
		};
		void internal_setNewModelParameters(Utility::InputFileReader& input, const Models::ModelParameters& modelParameters)
		{
			this->model->setNewModelParameters(modelParameters, Model::SYSTEM_TYPE);
			this->initialize_model(input, modelParameters);
		};

	public:
		virtual ~DetailModelConstructor() = default;

		DetailModelConstructor(Utility::InputFileReader& input, const Models::ModelParameters& modelParameters)
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