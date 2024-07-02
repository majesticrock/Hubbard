#include "ModeHandler.hpp"
#include <vector>
#include <filesystem>
#include "../Hubbard/Helper/SquareGeneral.hpp"
#include "../Hubbard/Helper/SquareXP.hpp"
#include "../Hubbard/Helper/DOSGeneral.hpp"
#include "../Hubbard/Helper/DOS_XP.hpp"
#include "../Hubbard/DensityOfStates/Square.hpp"
#include "../Hubbard/DensityOfStates/SimpleCubic.hpp"
#include <Utility/OutputConvenience.hpp>

using data_vector = std::vector<Hubbard::global_floating_type>;
const std::string BASE_FOLDER = "../../data/modes/";

std::unique_ptr<Hubbard::Helper::ModeHelper> ModeHandler::getHelper(Utility::InputFileReader& input, Hubbard::Models::ModelParameters& modelParameters) const
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

std::unique_ptr<Hubbard::Helper::ModeHelper> ModeHandler::getHelper(Utility::InputFileReader& input) const
{
	std::vector<double> model_params = input.getDoubleList("model_parameters");
	Hubbard::Models::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
		0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));
	return getHelper(input, modelParameters);
}

std::vector<std::string> ModeHandler::getFileComments(Utility::InputFileReader& input, Hubbard::Helper::ModeHelper* modeHelper) const
{
	using std:: to_string;

	std::vector<std::string> comments;
	comments.push_back("Used DOS: " + input.getString("use_DOS"));
	comments.push_back("Discretization: " + input.getString("k_discretization"));
	comments.push_back("Lattice type: " + input.getString("lattice_type"));
	comments.push_back("Total Gap: " + to_string(modeHelper->getModel().getTotalGapValue()));
	return comments;
}

void ModeHandler::execute(Utility::InputFileReader& input) const
{
	using std::to_string;

	data_vector oneParticleEnergies;
	std::vector<Hubbard::ResolventReturnData> resolvents;

	std::unique_ptr<Hubbard::Helper::ModeHelper> modeHelper{ getHelper(input) };

	modeHelper->getModel().getAllEnergies(oneParticleEnergies);
	resolvents = modeHelper->computeCollectiveModes();

	if (rank == 0) {
		const std::string output_folder{ getOutputFolder(input) + modelParameters.getFolderName() };
		std::cout << "Saving data to folder " << BASE_FOLDER + output_folder << std::endl;
		std::filesystem::create_directories(BASE_FOLDER + output_folder);
		const std::vector<std::string> comments = getFileComments(input, modeHelper.get());

		if (!resolvents.empty()) {
			for (size_t i = 0U; i < resolvents.size(); ++i)
			{
				resolvents[i].writeDataToFile(BASE_FOLDER + output_folder + "resolvent_" + resolvents[i].name, comments);
			}
		}
		else {
			std::cout << "Resolvent returned an empty vector." << std::endl;
		}

		Utility::saveData(oneParticleEnergies, BASE_FOLDER + output_folder + "one_particle.dat.gz", comments);
	}
}