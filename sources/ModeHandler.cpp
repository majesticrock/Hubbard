#include "ModeHandler.hpp"
#include <vector>
#include <filesystem>
#include "Hubbard/Helper/SquareGeneral.hpp"
#include "Hubbard/Helper/SquareXP.hpp"
#include "Hubbard/Helper/DOSGeneral.hpp"
#include "Hubbard/Helper/DOS_XP.hpp"
#include "Utility/OutputConvenience.hpp"
#include "Hubbard/DOSModels/BroydenDOS.hpp"
#include "Hubbard/DensityOfStates/Square.hpp"
#include "Hubbard/DensityOfStates/SimpleCubic.hpp"

using data_vector = std::vector<Hubbard::global_floating_type>;
const std::string BASE_FOLDER = "../../data/modes/";

std::unique_ptr<Hubbard::Helper::ModeHelper> ModeHandler::getHelper(Utility::InputFileReader& input, Hubbard::ModelParameters& modelParameters) const
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
	Hubbard::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
		0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));
	return getHelper(input, modelParameters);
}

void ModeHandler::execute(Utility::InputFileReader& input) const
{
	using std::to_string;

	std::vector<data_vector> reciever;
	data_vector oneParticleEnergies;
	Hubbard::global_floating_type totalGapValue;
	std::vector<Hubbard::ResolventReturnData> resolvents;

	std::unique_ptr<Hubbard::Helper::ModeHelper> modeHelper{ getHelper(input) };

	totalGapValue = modeHelper->getModel().getTotalGapValue();
	modeHelper->getModel().getAllEnergies(oneParticleEnergies);
	resolvents = modeHelper->computeCollectiveModes(reciever);

	if (rank == 0) {
		std::string output_folder{ getOutputFolder(input) + modelParameters.getFolderName() };
		std::cout << "Saving data to folder " << BASE_FOLDER + output_folder << std::endl;
		std::filesystem::create_directories(BASE_FOLDER + output_folder);

		std::vector<std::string> comments;
		comments.push_back("Used DOS: " + input.getString("use_DOS"));
		comments.push_back("Discretization: " + input.getString("k_discretization"));
		comments.push_back("Lattice type: " + input.getString("lattice_type"));
		comments.push_back("Total Gap: " + to_string(totalGapValue));

		if (!(reciever.empty())) {
			Utility::saveData(reciever, BASE_FOLDER + output_folder + ".dat.gz", comments);
		}
		if (resolvents.size() > 0U) {
			std::vector<std::string> names;
			if (input.getInt("start_basis_at") == -1) {
				names = { "phase_SC", "phase_CDW", "phase_AFM", "higgs_SC", "higgs_CDW", "higgs_AFM" };
			}
			else {
				names = { "higgs_SC_a", "higgs_SC_a+b", "higgs_SC_a+ib",
					"phase_SC_a", "phase_SC_a+b", "phase_SC_a+ib",
					"CDW_a", "CDW_a+b", "CDW_a+ib",
					"AFM_a", "AFM_a+b", "AFM_a+ib" };
			}

			for (size_t i = 0U; i < resolvents.size(); ++i)
			{
				resolvents[i].writeDataToFile(BASE_FOLDER + output_folder + "resolvent_" + names[i], comments);
			}
		}
		else {
			std::cout << "Resolvent returned an empty vector." << std::endl;
		}

		Utility::saveData(oneParticleEnergies, BASE_FOLDER + output_folder + "one_particle.dat.gz", comments);
	}
}