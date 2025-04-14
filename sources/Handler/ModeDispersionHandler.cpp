#include "ModeDispersionHandler.hpp"

#include <vector>
#include <array>
#include <filesystem>
#include <mrock/utility/info_to_json.hpp>
#include <mrock/info.h>
#include <nlohmann/json.hpp>
// File is generated on build by cmake
#include "../../build_header/info.h"

#include "../Hubbard/Helper/SquareGeneral.hpp"

const std::string BASE_FOLDER = "../../data/hubbard/";

void ModeDispersionHandler::execute(mrock::utility::InputFileReader& input) const
{
	using std::to_string;

	std::vector<double> model_params = input.getDoubleList("model_parameters");
	Hubbard::Models::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
		0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));

	std::vector<Hubbard::ResolventReturnData> resolvents;
	Hubbard::Helper::SquareGeneral modeHelper(input, modelParameters);

	if(eval_index < 0) {
		for (int i = 0; i < Hubbard::Constants::K_DISCRETIZATION; ++i)
		{
			mrock::utility::Numerics::join_data_wrapper(resolvents, modeHelper.compute_collective_modes());
			modeHelper.mode_momentum.x() += 1;
		}
		for (int i = 0; i < Hubbard::Constants::K_DISCRETIZATION; ++i)
		{
			mrock::utility::Numerics::join_data_wrapper(resolvents, modeHelper.compute_collective_modes());
			modeHelper.mode_momentum.y() += 1;
		}
		for (int i = 0; i < Hubbard::Constants::K_DISCRETIZATION; ++i)
		{
			mrock::utility::Numerics::join_data_wrapper(resolvents, modeHelper.compute_collective_modes());
			modeHelper.mode_momentum.x() -= 1;
			modeHelper.mode_momentum.y() -= 1;
		}
	} 
	else {
		modeHelper.mode_momentum = eval_point(eval_index);
		mrock::utility::Numerics::join_data_wrapper(resolvents, modeHelper.compute_collective_modes());
	}

	const std::string output_folder{ getOutputFolder(input) + modelParameters.getFolderName() };
	std::cout << "Saving data to folder " << BASE_FOLDER + output_folder << std::endl;
	std::filesystem::create_directories(BASE_FOLDER + output_folder);
	const std::vector<std::string> comments = getFileComments(input, &modeHelper);
	if (!resolvents.empty()) {
		nlohmann::json jResolvents = {
			{ "resolvents", resolvents },
			{ "time", mrock::utility::time_stamp() },
			{ "used_dos", input.getBool("use_DOS") },
			{ "discretization", input.getInt("k_discretization") },
			{ "lattice_type", input.getString("lattice_type") },
			{ "gap_parameters", modeHelper.getModel().getAttributes().selfconsistency_values },
			{ "total_gap", modeHelper.getModel().getTotalGapValue() },
			{ "continuum_boundaries", modeHelper.getModel().continuum_boundaries() },
			{ "T", modeHelper.getModel().temperature },
			{ "U", modeHelper.getModel().U },
			{ "V", modeHelper.getModel().V },
			{ "XP_basis", (input.getInt("start_basis_at") < 0 ? 1 : 0) },
			{ "start_ratio_cdw_sc", input.getDouble("ratio_CDW_SC") }
		};
		mrock::utility::saveString(jResolvents.dump(4), BASE_FOLDER + output_folder + std::to_string(eval_index) + "dispersions.json.gz");
	}
	else {
		std::cout << "Resolvent returned an empty vector." << std::endl;
	}

	// Generate metadata
	nlohmann::json info_json = mrock::utility::generate_json<Hubbard::info>("hubbard_");
	info_json.update(mrock::utility::generate_json<mrock::info>("mrock_"));
	mrock::utility::saveString(info_json.dump(4), BASE_FOLDER + output_folder + "metadata.json.gz");
}