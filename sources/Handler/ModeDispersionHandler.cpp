#include "ModeDispersionHandler.hpp"
#include "../Hubbard/Helper/SquareGeneral.hpp"
#include <vector>
#include <array>
#include <filesystem>

const std::string BASE_FOLDER = "../../data/hubbard/";

void ModeDispersionHandler::execute(Utility::InputFileReader& input) const
{
	using std::to_string;

	std::vector<double> model_params = input.getDoubleList("model_parameters");
	Hubbard::Models::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
		0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));

	std::vector<Hubbard::ResolventReturnData> resolvents;
	Hubbard::Helper::SquareGeneral modeHelper(input, modelParameters);

	const int TOTAL_EVAL_POINTS = 3 * Hubbard::Constants::K_DISCRETIZATION;
	const int EVAL_POINTS_PER_RANK = TOTAL_EVAL_POINTS / numberOfRanks;
	for (int i = 0; i < EVAL_POINTS_PER_RANK; ++i) {
		modeHelper.mode_momentum = eval_point(i + rank * EVAL_POINTS_PER_RANK);
		Utility::Numerics::join_data_wrapper(resolvents, modeHelper.computeCollectiveModes());
	}

	/*for (int i = 0; i < Hubbard::Constants::K_DISCRETIZATION; ++i)
	{
		Utility::Numerics::join_data_wrapper(resolvents, modeHelper.computeCollectiveModes());
		modeHelper.mode_momentum.x() += 1;
	}
	for (int i = 0; i < Hubbard::Constants::K_DISCRETIZATION; ++i)
	{
		Utility::Numerics::join_data_wrapper(resolvents, modeHelper.computeCollectiveModes());
		modeHelper.mode_momentum.y() += 1;
	}
	for (int i = 0; i < Hubbard::Constants::K_DISCRETIZATION; ++i)
	{
		Utility::Numerics::join_data_wrapper(resolvents, modeHelper.computeCollectiveModes());
		modeHelper.mode_momentum.x() -= 1;
		modeHelper.mode_momentum.y() -= 1;
	}*/

	const std::string output_folder{ getOutputFolder(input) + modelParameters.getFolderName() };
	std::cout << "Saving data to folder " << BASE_FOLDER + output_folder << std::endl;
	std::filesystem::create_directories(BASE_FOLDER + output_folder);
	const std::vector<std::string> comments = getFileComments(input, &modeHelper);
	if (!resolvents.empty()) {
		nlohmann::json jResolvents = {
			{ "resolvents", resolvents },
			{ "time", Utility::time_stamp() },
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
		Utility::saveString(jResolvents.dump(4), BASE_FOLDER + output_folder + "dispersions.json.gz");
	}
	else {
		std::cout << "Resolvent returned an empty vector." << std::endl;
	}
}