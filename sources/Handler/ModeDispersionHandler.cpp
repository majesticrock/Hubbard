#include "ModeDispersionHandler.hpp"
#include "../Hubbard/Helper/SquareXP.hpp"
#include <vector>
#include <array>
#include <filesystem>

const std::string BASE_FOLDER = "../../data/dispersions/";

void ModeDispersionHandler::execute(Utility::InputFileReader& input) const
{
    using std::to_string;

    std::vector<double> model_params = input.getDoubleList("model_parameters");
	Hubbard::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
		0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));

    std::vector<Hubbard::ResolventReturnData> resolvents;
    Hubbard::Helper::SquareXP modeHelper(input, modelParameters);

    const std::vector<std::string> names{ "phase_SC", "phase_CDW", "phase_AFM", "phase_AFM_trans", "higgs_SC", "higgs_CDW", "higgs_AFM", "higgs_AFM_trans" };
    const std::string output_folder{ getOutputFolder(input) + modelParameters.getFolderName() };
    std::filesystem::create_directories(BASE_FOLDER + output_folder);

    const int TOTAL_EVAL_POINTS = 3 * Hubbard::Constants::K_DISCRETIZATION;
    for(int i = 0; i < TOTAL_EVAL_POINTS; ++i)
    {
        modeHelper.computeCollectiveModes();
        if(i == TOTAL_EVAL_POINTS - 1){
            resolvents = modeHelper.computeCollectiveModes();
        }
    }

    if (!resolvents.empty()) {
        const std::vector<std::string> comments = getFileComments(input, &modeHelper);
        for (size_t i = 0U; i < resolvents.size(); ++i) {
			resolvents[i].writeDataToFile(BASE_FOLDER + output_folder + "resolvent_" + names[i], comments);
		}
	}
	else {
		std::cout << "Resolvent returned an empty vector." << std::endl;
	}
}