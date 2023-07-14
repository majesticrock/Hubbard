#ifdef _DEBUG
#define _NO_MPI
#endif

#include <omp.h>
#ifndef _NO_MPI
#include <mpi.h>
#endif

#include <filesystem>
#include "Hubbard/Helper/GeneralBasis.hpp"
#include "Hubbard/Helper/PhaseHelper.hpp"
#include "Hubbard/Helper/XPModes.hpp"
#include "Hubbard/SquareLattice/SquareTripletPairing.hpp"
#include "Hubbard/DOSModels/BroydenDOS.hpp"
#include "Hubbard/DensityOfStates/Square.hpp"

using Hubbard::Helper::data_vector;

int Hubbard::Constants::K_DISCRETIZATION = 100;
int Hubbard::Constants::BASIS_SIZE = 10000;
constexpr int NUMBER_OF_PARAMETERS = 8;
constexpr int NUMBER_OF_GAP_VALUES = NUMBER_OF_PARAMETERS - 2; // The Fock parameters are not important
std::vector<std::string> Hubbard::Constants::option_list = { "T", "U", "V" };

std::ostream& operator<<(std::ostream& os, const Hubbard::ModelParameters& mp) {
	os << mp.temperature << "\t" << mp.U << "\t" << mp.V;
	return os;
}

int main(int argc, char** argv)
{
#ifndef _NO_MPI
	if (argc < 2) {
		std::cerr << "Invalid number of arguments: Use mpirun -n <threads> <path_to_executable> <configfile>" << std::endl;
		return -1;
	}
	// First call MPI_Init
	MPI_Init(&argc, &argv);

	// Get my rank and the number of ranks
	int rank, numberOfRanks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfRanks);

#else
	int rank = 0;
	int numberOfRanks = 1;
#endif
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	if (rank == 0) {
		std::cout << "Using parameter file " << argv[1] << std::endl;
	}
	Utility::InputFileReader input(argv[1]);
	Hubbard::Constants::K_DISCRETIZATION = input.getInt("k_discretization");
	Hubbard::Constants::BASIS_SIZE = 4 * Hubbard::Constants::K_DISCRETIZATION * Hubbard::Constants::K_DISCRETIZATION;
	// Setup the parameters T, U, V
	std::vector<double> model_params = input.getDoubleList("model_parameters");

	if (input.getString("compute_what") == "test") {
		std::chrono::steady_clock::time_point test_b = std::chrono::steady_clock::now();
		std::chrono::steady_clock::time_point test_e;
		Hubbard::ModelParameters mP(model_params[0], model_params[1], model_params[2], 0, 0, "", "");
		
		//------------------------------------------------------------//

		Hubbard::DOSModels::BroydenDOS<Hubbard::DensityOfStates::Square> model3(mP);
		test_b = std::chrono::steady_clock::now();
		model3.computePhases({ true, true }).print();
		std::cout << "Free energy = " << model3.freeEnergyPerSite() << std::endl;

		test_e = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
		return 0;


		Hubbard::SquareLattice::HubbardCDW model(mP);
		//Hubbard::Constants::BASIS_SIZE = 2 * Hubbard::Constants::K_DISCRETIZATION;
		//Hubbard::ChainLattice::TripletPairingIterative model(mP);
		//model.computePhases({false, true}).print();
		std::cout << "Free energy = " << model.freeEnergyPerSite() << std::endl;

		test_e = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
		std::cout << "\n\n" << std::endl;

#ifndef _NO_MPI
		//return MPI_Finalize();
#else
		//return 0;
#endif // !_NO_MPI

		Hubbard::SquareLattice::UsingBroyden model2(mP);
		test_b = std::chrono::steady_clock::now();
		//model2.computePhases({false, true}).print();
		std::cout << "Free energy = " << model2.freeEnergyPerSite() << std::endl;

		test_e = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
	}
	else if (input.getString("compute_what") == "phases") {
		// Setup the number of steps
		int GLOBAL_IT_STEPS = input.getInt("global_iterator_steps");
		int FIRST_IT_STEPS = GLOBAL_IT_STEPS / numberOfRanks;
		double GLOBAL_IT_LIMS[2] = { 0, input.getDouble("global_iterator_upper_limit") };
		std::vector<std::string> option_list = { "T", "U", "V" };
		double FIRST_IT_RANGE = 0;
		double FIRST_IT_MIN = 0;
		for (int i = 0; i < option_list.size(); i++)
		{
			if (input.getString("global_iterator_type") == option_list[i]) {
				GLOBAL_IT_LIMS[0] = model_params[i];
				FIRST_IT_RANGE = (GLOBAL_IT_LIMS[1] - GLOBAL_IT_LIMS[0]) / numberOfRanks;
				FIRST_IT_MIN = GLOBAL_IT_LIMS[0] + rank * FIRST_IT_RANGE;
				model_params[i] = FIRST_IT_MIN;
			}
		}

		int SECOND_IT_STEPS = input.getInt("second_iterator_steps");
		double SECOND_IT_MIN = 0, SECOND_IT_MAX = input.getDouble("second_iterator_upper_limit");

		for (int i = 0; i < option_list.size(); i++)
		{
			if (input.getString("second_iterator_type") == option_list[i]) {
				SECOND_IT_MIN = model_params[i];
			}
		}

		std::vector<data_vector> local_data(NUMBER_OF_PARAMETERS, data_vector(FIRST_IT_STEPS * SECOND_IT_STEPS));

		Hubbard::Helper::PhaseHelper phaseHelper(input, rank, numberOfRanks);
		phaseHelper.compute_crude(local_data);

		std::vector<data_vector> recieve_data(NUMBER_OF_PARAMETERS, data_vector(GLOBAL_IT_STEPS * SECOND_IT_STEPS));

#ifndef _NO_MPI
		for (size_t i = 0; i < NUMBER_OF_PARAMETERS; i++)
		{
			MPI_Allgather(local_data[i].data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE,
				recieve_data[i].data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, MPI_COMM_WORLD);
		}
#endif
		if (rank == 0) {
			std::vector<std::string> comments;
			comments.push_back(input.getString("second_iterator_type") + "_MIN=" + std::to_string(SECOND_IT_MIN)
				+ "   " + input.getString("second_iterator_type") + "_MAX=" + std::to_string(SECOND_IT_MAX));

			comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
				+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

			std::string output_folder = input.getString("output_folder");
			std::filesystem::create_directories("../../data/phases/" + output_folder);

			Utility::saveData(recieve_data[0], SECOND_IT_STEPS, "../../data/phases/" + output_folder + "cdw.dat.gz", comments);
			Utility::saveData(recieve_data[1], SECOND_IT_STEPS, "../../data/phases/" + output_folder + "afm.dat.gz", comments);
			Utility::saveData(recieve_data[2], SECOND_IT_STEPS, "../../data/phases/" + output_folder + "sc.dat.gz", comments);
			Utility::saveData(recieve_data[3], SECOND_IT_STEPS, "../../data/phases/" + output_folder + "gamma_sc.dat.gz", comments);
			Utility::saveData(recieve_data[4], SECOND_IT_STEPS, "../../data/phases/" + output_folder + "xi_sc.dat.gz", comments);
			Utility::saveData(recieve_data[5], SECOND_IT_STEPS, "../../data/phases/" + output_folder + "eta.dat.gz", comments);

			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << "Crude computation time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		}

		if (input.getBool("improved_boundaries")) {
			std::vector<data_vector> local(NUMBER_OF_GAP_VALUES);
			// We stick to int (and ignore size_t vs. int mismatch) due to MPI restrictions
			std::vector<int> sizes(NUMBER_OF_GAP_VALUES);
			std::vector<std::vector<int>> all_sizes(NUMBER_OF_GAP_VALUES);
			std::vector<int> totalSizes(NUMBER_OF_GAP_VALUES, 0);

			for (size_t i = 0U; i < NUMBER_OF_GAP_VALUES; ++i)
			{
				phaseHelper.findSingleBoundary(recieve_data, local[i], i);
				sizes[i] = local[i].size();
			}
			std::vector < data_vector> recieve_boundaries(NUMBER_OF_GAP_VALUES);
			if (rank == 0) {
				for (size_t i = 0U; i < NUMBER_OF_GAP_VALUES; ++i)
				{
					all_sizes[i].resize(numberOfRanks);
				}
			}
#ifndef _NO_MPI
			for (size_t i = 0; i < NUMBER_OF_GAP_VALUES; i++)
			{
				MPI_Gather(&(sizes[i]), 1, MPI_INT, all_sizes[i].data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
				std::vector<int> displacements(numberOfRanks, 0);

				if (rank == 0) {
					for (const auto& s : all_sizes[i])
					{
						totalSizes[i] += s;
					}
					for (size_t j = 1; j < numberOfRanks; j++)
					{
						displacements[j] = displacements[j - 1] + all_sizes[i][j - 1];
					}
				}

				recieve_boundaries[i].resize(totalSizes[i]);
				MPI_Gatherv(local[i].data(), sizes[i], MPI_DOUBLE, recieve_boundaries[i].data(), all_sizes[i].data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
#endif
			if (rank == 0) {
				std::vector<std::string> comments;
				comments.push_back(input.getString("second_iterator_type") + "_MIN=" + std::to_string(SECOND_IT_MIN)
					+ "   " + input.getString("second_iterator_type") + "_MAX=" + std::to_string(SECOND_IT_MAX));

				comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
					+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

				std::string output_folder = input.getString("output_folder");
				std::filesystem::create_directories("../../data/phases/" + output_folder);

				std::string names[] = { "cdw", "afm", "sc", "gamma_sc", "xi_sc", "eta" };
				for (size_t i = 0U; i < NUMBER_OF_GAP_VALUES; ++i)
				{
					const size_t n = recieve_boundaries[i].size() / 2;
					std::vector<std::vector<double>> buffer(2, std::vector<double>(n));
					for (size_t j = 0U; j < recieve_boundaries[i].size(); j += 2)
					{
						buffer[0][j / 2] = recieve_boundaries[i][j];
						buffer[1][j / 2] = recieve_boundaries[i][j + 1];
					}

					Utility::saveData(buffer, "../../data/phases/" + output_folder + "boundaries_" + names[i] + ".dat.gz", comments);
				}
			}
		}
	}
	else if (input.getString("compute_what") == "modes") {
		omp_set_num_threads(8);

		std::vector<data_vector> reciever;
		std::vector<data_vector> oneParticleEnergies;
		double totalGapValue;
		std::vector<Hubbard::Resolvent_L> resolvents;

		std::unique_ptr<Hubbard::Helper::ModeHelper> modeHelper;
		if (input.getInt("start_basis_at") == -1) {
			modeHelper = std::make_unique<Hubbard::Helper::XPModes>(input);
		}
		else {
			modeHelper = std::make_unique<Hubbard::Helper::GeneralBasis>(input);
		}

		totalGapValue = modeHelper->getModel().getTotalGapValue();
		modeHelper->getModel().getAllEnergies(oneParticleEnergies);
		resolvents = modeHelper->computeCollectiveModes(reciever);

		if (rank == 0) {
			Hubbard::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
				0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));
			std::string output_folder = input.getString("output_folder") + modelParameters.getFileName();
			std::filesystem::create_directories("../../data/" + output_folder);

			std::vector<std::string> comments;
			comments.push_back("Total Gap=" + std::to_string(totalGapValue));
			if (!(reciever.empty())) {
				Utility::saveData(reciever, "../../data/" + output_folder + ".dat.gz", comments);
			}
			if (resolvents.size() > 0) {
				std::vector<std::string> names;
				if (input.getInt("start_basis_at") == -1) {
					names = { "phase_SC", "phase_CDW", "phase_AFM", "higgs_SC", "higgs_CDW", "higgs_AFM" };
				}
				else {
					names = { "higgs_sc_a", "higgs_sc_a+b", "higgs_sc_a+ib",
						"phase_sc_a", "phase_sc_a+b", "phase_sc_a+ib",
						"cdw_a", "cdw_a+b", "cdw_a+ib",
						"afm_a", "afm_a+b", "afm_a+ib" };
				}

				for (size_t i = 0; i < resolvents.size(); i++)
				{
					resolvents[i].writeDataToFile("../../data/" + output_folder + "resolvent_" + names[i]);
				}
			}
			else {
				std::cout << "Resolvent returned an empty vector." << std::endl;
			}
			comments.pop_back();
			Utility::saveData(oneParticleEnergies, "../../data/" + output_folder + "one_particle.dat.gz", comments);
		}
	}

	if (rank == 0) {
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	}
#ifndef _NO_MPI
	return MPI_Finalize();
#else
	return 0;
#endif
}