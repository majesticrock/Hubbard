#ifdef _DEBUG
#define _NO_MPI
#endif

#include <omp.h>
#ifndef _NO_MPI
#include <mpi.h>
#endif

#include "PhaseHandler.hpp"
#include "Hubbard/Constants.hpp"
#include "Hubbard/Helper/PhaseHelper.hpp"
#include "Utility/OutputConvenience.hpp"
#include <vector>
#include <chrono>
#include <filesystem>

using data_vector = std::vector<Hubbard::global_floating_type>;
constexpr int NUMBER_OF_PARAMETERS = 8;
constexpr int NUMBER_OF_GAP_VALUES = NUMBER_OF_PARAMETERS - 2; // The Fock parameters are not important
const std::string BASE_FOLDER = "../../data/phases/";

void PhaseHandler::execute(Utility::InputFileReader& input) const
{
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	std::vector<double> model_params{ input.getDoubleList("model_parameters") };
	// Setup the number of steps
	int GLOBAL_IT_STEPS = input.getInt("global_iterator_steps");
	int FIRST_IT_STEPS = GLOBAL_IT_STEPS / numberOfRanks;
	double GLOBAL_IT_LIMS[2] = { 0, input.getDouble("global_iterator_upper_limit") };
	double FIRST_IT_RANGE = 0;
	double FIRST_IT_MIN = 0;
	for (int i = 0; i < Hubbard::Constants::option_list.size(); i++)
	{
		if (input.getString("global_iterator_type") == Hubbard::Constants::option_list[i]) {
			GLOBAL_IT_LIMS[0] = model_params[i];
			FIRST_IT_RANGE = (GLOBAL_IT_LIMS[1] - GLOBAL_IT_LIMS[0]) / numberOfRanks;
			FIRST_IT_MIN = GLOBAL_IT_LIMS[0] + rank * FIRST_IT_RANGE;
			model_params[i] = FIRST_IT_MIN;
		}
	}

	int SECOND_IT_STEPS = input.getInt("second_iterator_steps");
	double SECOND_IT_MIN = 0, SECOND_IT_MAX = input.getDouble("second_iterator_upper_limit");

	for (size_t i = 0U; i < Hubbard::Constants::option_list.size(); ++i)
	{
		if (input.getString("second_iterator_type") == Hubbard::Constants::option_list[i]) {
			SECOND_IT_MIN = model_params[i];
		}
	}

	std::vector<data_vector> local_data(NUMBER_OF_PARAMETERS, data_vector(FIRST_IT_STEPS * SECOND_IT_STEPS));
	Hubbard::Helper::PhaseHelper phaseHelper(input, rank, numberOfRanks);
	phaseHelper.compute_crude(local_data);

	std::vector<data_vector> local_coexitence_data(3, data_vector(FIRST_IT_STEPS));
	phaseHelper.coexistence_AFM_CDW(local_coexitence_data);

	std::vector<data_vector> recieve_data(NUMBER_OF_PARAMETERS, data_vector(GLOBAL_IT_STEPS * SECOND_IT_STEPS));
	std::vector<data_vector> recieve_coexitence_data(3, data_vector(GLOBAL_IT_STEPS));

#ifndef _NO_MPI
	for (size_t i = 0U; i < NUMBER_OF_PARAMETERS; ++i)
	{
		MPI_Allgather(local_data[i].data(), FIRST_IT_STEPS * SECOND_IT_STEPS, _MPI_RETURN_TYPE,
			recieve_data[i].data(), FIRST_IT_STEPS * SECOND_IT_STEPS, _MPI_RETURN_TYPE, MPI_COMM_WORLD);
	}
	for (size_t i = 0U; i < 3; ++i){
		MPI_Gather(local_coexitence_data[i].data(), FIRST_IT_STEPS, MPI_DOUBLE, 
			recieve_coexitence_data[i].data(), FIRST_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
#else
	for (size_t i = 0U; i < NUMBER_OF_PARAMETERS; ++i)
	{
		recieve_data[i] = local_data[i];
	}
	for (size_t i = 0U; i < 3; ++i){
		recieve_coexitence_data[i] = local_coexitence_data[i];
	}
#endif
	std::string output_folder{ getOutputFolder(input) };
	if (rank == 0) {
		std::vector<std::string> comments;
		comments.push_back(input.getString("second_iterator_type") + "_MIN=" + std::to_string(SECOND_IT_MIN)
			+ "   " + input.getString("second_iterator_type") + "_MAX=" + std::to_string(SECOND_IT_MAX));

		comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
			+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

		comments.push_back("Used DOS: " + input.getString("use_DOS"));
		comments.push_back("Discretization: " + input.getString("k_discretization"));
		comments.push_back("Lattice type: " + input.getString("lattice_type"));

		std::cout << "Saving data to folder " << BASE_FOLDER + output_folder << std::endl;
		std::filesystem::create_directories(BASE_FOLDER + output_folder);

		Utility::saveData(recieve_data[0], SECOND_IT_STEPS, BASE_FOLDER + output_folder + "cdw.dat.gz", comments);
		Utility::saveData(recieve_data[1], SECOND_IT_STEPS, BASE_FOLDER + output_folder + "afm.dat.gz", comments);
		Utility::saveData(recieve_data[2], SECOND_IT_STEPS, BASE_FOLDER + output_folder + "sc.dat.gz", comments);
		Utility::saveData(recieve_data[3], SECOND_IT_STEPS, BASE_FOLDER + output_folder + "gamma_sc.dat.gz", comments);
		Utility::saveData(recieve_data[4], SECOND_IT_STEPS, BASE_FOLDER + output_folder + "xi_sc.dat.gz", comments);
		Utility::saveData(recieve_data[5], SECOND_IT_STEPS, BASE_FOLDER + output_folder + "eta.dat.gz", comments);

		Utility::saveData(recieve_coexitence_data, BASE_FOLDER + output_folder + "coexistence_afm_cdw.dat.gz", comments);

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
		std::vector<data_vector> recieve_boundaries(NUMBER_OF_GAP_VALUES);
		if (rank == 0) {
			for (size_t i = 0U; i < NUMBER_OF_GAP_VALUES; ++i)
			{
				all_sizes[i].resize(numberOfRanks);
			}
		}
#ifndef _NO_MPI
		for (size_t i = 0U; i < NUMBER_OF_GAP_VALUES; ++i)
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
			MPI_Gatherv(local[i].data(), sizes[i], _MPI_RETURN_TYPE, recieve_boundaries[i].data(), all_sizes[i].data(), displacements.data(), _MPI_RETURN_TYPE, 0, MPI_COMM_WORLD);
		}
#else
		for (size_t i = 0U; i < NUMBER_OF_GAP_VALUES; ++i)
		{
			recieve_boundaries[i] = local[i];
		}
#endif
		if (rank == 0) {
			std::vector<std::string> comments;
			comments.push_back(input.getString("second_iterator_type") + "_MIN=" + std::to_string(SECOND_IT_MIN)
				+ "   " + input.getString("second_iterator_type") + "_MAX=" + std::to_string(SECOND_IT_MAX));

			comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
				+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

			std::filesystem::create_directories(BASE_FOLDER + output_folder);

			std::string names[] = { "cdw", "afm", "sc", "gamma_sc", "xi_sc", "eta" };
			for (size_t i = 0U; i < NUMBER_OF_GAP_VALUES; ++i)
			{
				const size_t n = recieve_boundaries[i].size() / 2;
				std::vector<data_vector> buffer(2, data_vector(n));
				for (size_t j = 0U; j < recieve_boundaries[i].size(); j += 2)
				{
					buffer[0][j / 2] = recieve_boundaries[i][j];
					buffer[1][j / 2] = recieve_boundaries[i][j + 1];
				}

				Utility::saveData(buffer, BASE_FOLDER + output_folder + "boundaries_" + names[i] + ".dat.gz", comments);
			}
		}
		}
	}