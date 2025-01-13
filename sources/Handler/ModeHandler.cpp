#include "ModeHandler.hpp"

#include <vector>
#include <filesystem>

#include <mrock/utility/OutputConvenience.hpp>
#include <mrock/utility/info_to_json.hpp>
#include <mrock/info.h>
#include <nlohmann/json.hpp>
// File is generated on build by cmake
#include "../../build/info.h"

#include "../Hubbard/Helper/SquareGeneral.hpp"
#include "../Hubbard/Helper/SquareXP.hpp"
#include "../Hubbard/Helper/DOSGeneral.hpp"
#include "../Hubbard/Helper/DOS_XP.hpp"
#include "../Hubbard/DensityOfStates/Square.hpp"
#include "../Hubbard/DensityOfStates/SimpleCubic.hpp"


using data_vector = std::vector<Hubbard::global_floating_type>;
const std::string BASE_FOLDER = "../../data/hubbard/";

std::unique_ptr<Hubbard::Helper::ModeHelper> ModeHandler::getHelper(mrock::utility::InputFileReader& input, Hubbard::Models::ModelParameters& modelParameters) const
{
	if (input.getInt("start_basis_at") == -1) {
		if (input.getBool("use_DOS")) {
			if (input.getString("lattice_type") == "square") {
				return std::make_unique<Hubbard::Helper::DOS_XP<Hubbard::DensityOfStates::Square>>(input, modelParameters);
			}
			else if (input.getString("lattice_type") == "cube") {
				return std::make_unique<Hubbard::Helper::DOS_XP<Hubbard::DensityOfStates::SimpleCubic>>(input, modelParameters);
			}
			else {
				throw std::runtime_error("Could not find lattice_type: " + input.getString("lattice_type"));
			}
		}
		else {
			return std::make_unique<Hubbard::Helper::SquareXP>(input, modelParameters);
		}
	}
	else {
		if (input.getBool("use_DOS")) {
			if (input.getString("lattice_type") == "square") {
				return std::make_unique<Hubbard::Helper::DOSGeneral<Hubbard::DensityOfStates::Square>>(input, modelParameters);
			}
			else if (input.getString("lattice_type") == "cube") {
				return std::make_unique<Hubbard::Helper::DOSGeneral<Hubbard::DensityOfStates::SimpleCubic>>(input, modelParameters);
			}
			else {
				throw std::runtime_error("Could not find lattice_type: " + input.getString("lattice_type"));
			}
		}
		else {
			return std::make_unique<Hubbard::Helper::SquareGeneral>(input, modelParameters);
		}
	}
	return nullptr;
}

std::unique_ptr<Hubbard::Helper::ModeHelper> ModeHandler::getHelper(mrock::utility::InputFileReader& input) const
{
	std::vector<double> model_params = input.getDoubleList("model_parameters");
	Hubbard::Models::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
		0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));
	return getHelper(input, modelParameters);
}

std::vector<std::string> ModeHandler::getFileComments(mrock::utility::InputFileReader& input, Hubbard::Helper::ModeHelper* modeHelper) const
{
	using std::to_string;

	std::vector<std::string> comments;
	comments.push_back("Used DOS: " + input.getString("use_DOS"));
	comments.push_back("Discretization: " + input.getString("k_discretization"));
	comments.push_back("Lattice type: " + input.getString("lattice_type"));
	comments.push_back("Total Gap: " + to_string(modeHelper->getModel().getTotalGapValue()));
	return comments;
}

void ModeHandler::execute(mrock::utility::InputFileReader& input) const
{
	using std::to_string;

	data_vector oneParticleEnergies;
	std::unique_ptr<Hubbard::Helper::ModeHelper> modeHelper{ getHelper(input) };
	modeHelper->getModel().getAllEnergies(oneParticleEnergies);
	std::vector<Hubbard::ResolventReturnData> resolvents = modeHelper->compute_collective_modes();

	if (rank == 0) {
		const std::string output_folder{ getOutputFolder(input) + modelParameters.getFolderName() };
		std::cout << "Saving data to folder " << BASE_FOLDER + output_folder << std::endl;
		std::filesystem::create_directories(BASE_FOLDER + output_folder);
		const std::vector<std::string> comments = getFileComments(input, modeHelper.get());

		if (!resolvents.empty()) {
			nlohmann::json jResolvents = {
				{ "resolvents", resolvents },
				{ "time", mrock::utility::time_stamp() },
				{ "used_dos", input.getBool("use_DOS") },
				{ "discretization", input.getInt("k_discretization") },
				{ "lattice_type", input.getString("lattice_type") },
				{ "total_gap", modeHelper->getModel().getTotalGapValue() },
				{ "continuum_boundaries", modeHelper->getModel().continuum_boundaries() },
				{ "T", modeHelper->getModel().temperature },
				{ "U", modeHelper->getModel().U },
				{ "V", modeHelper->getModel().V },
				{ "XP_basis", (input.getInt("start_basis_at") < 0 ? 1 : 0) },
				{ "start_ratio_cdw_sc", input.getDouble("ratio_CDW_SC") }
			};
			mrock::utility::saveString(jResolvents.dump(4), BASE_FOLDER + output_folder + "resolvents.json.gz");
		}
		else {
			std::cout << "Resolvent returned an empty vector." << std::endl;
		}

		// Generate metadata
		nlohmann::json info_json = mrock::utility::generate_json<Hubbard::info>("hubbard_");
		info_json.update(mrock::utility::generate_json<mrock::info>("mrock_"));
		mrock::utility::saveString(info_json.dump(4), BASE_FOLDER + output_folder + "metadata.json.gz");
	}
}