#include "ModeDispersionHandler.hpp"
#include "../Hubbard/Helper/SquareXP.hpp"
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
	Hubbard::Helper::SquareXP modeHelper(input, modelParameters);

	const int TOTAL_EVAL_POINTS = 3 * Hubbard::Constants::K_DISCRETIZATION;
	for(int i = 0; i < Hubbard::Constants::K_DISCRETIZATION; ++i)
	{
		modeHelper.computeCollectiveModes();
		modeHelper.mode_momentum.x() += 1;
	}
	for (int i = 0; i < Hubbard::Constants::K_DISCRETIZATION; ++i)
	{
		modeHelper.computeCollectiveModes();
		modeHelper.mode_momentum.y() += 1;
	}
	for (int i = 0; i < Hubbard::Constants::K_DISCRETIZATION; ++i)
	{
		modeHelper.computeCollectiveModes();
		modeHelper.mode_momentum.x() -= 1;
		modeHelper.mode_momentum.y() -= 1;
	}

	const std::string output_folder{ getOutputFolder(input) + modelParameters.getFolderName() };
	std::cout << "Saving data to folder " << BASE_FOLDER + output_folder << std::endl;
	std::filesystem::create_directories(BASE_FOLDER + output_folder);
	const std::vector<std::string> comments = getFileComments(input, &modeHelper);
	if (!resolvents.empty()) {
		nlohmann::json jResolvents = {
			{ "resolvents", resolvents },
			{ "time", Utility::time_stamp() },
			{ "Used DOS", input.getBool("use_DOS") },
			{ "Discretization", input.getInt("k_discretization") },
			{ "Lattice type", input.getString("lattice_type") },
			{ "Total Gap", modeHelper.getModel().getTotalGapValue() },
			{ "Continuum Boundaries", modeHelper.getModel().continuum_boundaries() },
			{ "T", modeHelper.getModel().temperature }, 
			{ "U", modeHelper.getModel().U }, 
			{ "V", modeHelper.getModel().V }
		};
		Utility::saveString(jResolvents.dump(4), BASE_FOLDER + output_folder + "resolvents.json.gz");
	}
	else {
		std::cout << "Resolvent returned an empty vector." << std::endl;
	}
}