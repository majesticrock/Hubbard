#pragma once
#include <memory>
#include "../BaseModel.hpp"
#include "../../Utility/InputFileReader.hpp"
#include "../SquareLattice/UsingBroyden.hpp"
#include <map>
#include "../../../../FermionCommute/sources/WickTerm.hpp"

namespace Hubbard::Helper {
	template <class Model>
	class DetailModelConstructor {
	protected:
		std::vector<ValueArray> expecs{};
		ValueArray sum_of_all;
		std::unique_ptr<Model> model{};

		const std::map<std::string, int> wick_map = { {"n", 0}, {"g", 1}, {"f", 2}, {"\\eta", 3} };
		const std::map<std::string, int> wick_spin_offset = { {"\\uparrow", 0}, {"\\downarrow", 4}, {"\\sigma", 6} };

		inline complex_prec getSumOfAll(const SymbolicOperators::WickOperator& op, int cos_modulation = 0) const {
			auto it = wick_map.find(op.type);
			if (it == wick_map.end()) throw std::invalid_argument("Term type not recognized: " + op.type);

			int index = it->second;
			if (op.type == "g" || op.type == "n") {
				auto jt = wick_spin_offset.find(op.indizes[0]);
				if (jt == wick_map.end()) throw std::runtime_error("Something went wrong while looking up the spin indizes.");
				index += jt->second;
			}

			if (op.isDaggered) return std::conj(sum_of_all(index, cos_modulation));
			return sum_of_all(index, cos_modulation);
		};

	public:
		virtual ~DetailModelConstructor() = default;
		DetailModelConstructor(Utility::InputFileReader& input)
		{
			std::vector<double> model_params = input.getDoubleList("model_parameters");
			Hubbard::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
				0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));

			model = std::make_unique<Model>(modelParameters);
			if (input.getString("ratio_CDW_SC") != "-1") {
				model->set_CDW_SC_ratio(input.getDouble("ratio_CDW_SC"));
			}
			ModelAttributes<global_floating_type> result = model->computePhases();

			if (model_params[1] > 0 && model_params[2] > 0) { // only for U>0 and V>0
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
					copy = model_copy_ptr->computePhases(NoWarning);
					if (copy.converged) {
						if (model_copy_ptr->freeEnergyPerSite() <= model->freeEnergyPerSite()) {
							model = std::move(model_copy_ptr);
						}
					}
				}
			}
			model->getAttributes().print();
			model->computeExpectationValues(expecs, sum_of_all);
		};

		DetailModelConstructor(std::unique_ptr<Model>&& model_ptr) : model{ std::move(model_ptr) }
		{
			model->computeExpectationValues(expecs, sum_of_all);
		};
	};
}