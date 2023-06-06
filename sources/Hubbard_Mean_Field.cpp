#define _USE_MATH_DEFINES

#ifdef _DEBUG
#define _NO_MPI
#endif

#include <omp.h>
#ifndef _NO_MPI
#include <mpi.h>
#endif

#include <string>
#include <iostream>
#include <filesystem>
#include <chrono>
#include <memory>
#include <cmath>
#include <limits>

#include "PhaseHelper.hpp"
#include "Utility/InputFileReader.hpp"
#include "Hubbard/BasicHubbardModel.hpp"
#include "Hubbard/HubbardCDW.hpp"
#include "Hubbard/UsingBroyden.hpp"
#include "Hubbard/Constants.hpp"
#include "Utility/OutputConvenience.hpp"

int Hubbard::Constants::K_DISCRETIZATION = 100;

std::ostream& operator<<(std::ostream& os, const Hubbard::Model::ModelParameters& mp) {
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
	// Setup the parameters T, U, V
	std::vector<double> model_params = input.getDoubleList("model_parameters");

	if (input.getString("compute_what") == "test") {
		Hubbard::Model::ModelParameters mP(model_params[0], model_params[1], model_params[2], 0, 0, "", "");
		Hubbard::HubbardCDW model(mP, 0, 0);

		std::chrono::steady_clock::time_point test_b = std::chrono::steady_clock::now();
		model.computePhases(true).print();
		std::chrono::steady_clock::time_point test_e = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;
		std::cout << "\n\n" << std::endl;
		return MPI_Finalize();
		Hubbard::UsingBroyden model2(mP, 0, 0);
		test_b = std::chrono::steady_clock::now();
		model2.computePhases(true).print();
		test_e = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;

		std::vector<std::vector<double>> energies;
		model2.getEnergies(energies, 1);
		Utility::saveData_boost(energies, "../../data/energies.dat.gz");
		return MPI_Finalize();
	}

	// Setup the number of steps
	int GLOBAL_IT_STEPS = input.getInt("global_iterator_steps");
	int FIRST_IT_STEPS = GLOBAL_IT_STEPS / numberOfRanks;
	double GLOBAL_IT_LIMS[2] = { 0, input.getDouble("global_iterator_upper_limit") };
	std::vector<std::string> option_list = { "T", "U", "V" };
	double FIRST_IT_RANGE = 0;
	double FIRST_IT_MIN = 0, FIRST_IT_MAX = 0;
	for (int i = 0; i < option_list.size(); i++)
	{
		if (input.getString("global_iterator_type") == option_list[i]) {
			GLOBAL_IT_LIMS[0] = model_params[i];
			FIRST_IT_RANGE = (GLOBAL_IT_LIMS[1] - GLOBAL_IT_LIMS[0]) / numberOfRanks;
			FIRST_IT_MIN = GLOBAL_IT_LIMS[0] + rank * FIRST_IT_RANGE;
			FIRST_IT_MAX = FIRST_IT_MIN + FIRST_IT_RANGE;
			model_params[i] = FIRST_IT_MIN;
		}
	}

	if (input.getString("compute_what") == "phases") {
		int SECOND_IT_STEPS = input.getInt("second_iterator_steps");
		double SECOND_IT_MIN = 0, SECOND_IT_MAX = input.getDouble("second_iterator_upper_limit");

		for (int i = 0; i < option_list.size(); i++)
		{
			if (input.getString("second_iterator_type") == option_list[i]) {
				SECOND_IT_MIN = model_params[i];
			}
		}

		Hubbard::Model::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
			(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, (SECOND_IT_MAX - SECOND_IT_MIN) / SECOND_IT_STEPS,
			input.getString("global_iterator_type"), input.getString("second_iterator_type"));

		data_vector data_cdw(FIRST_IT_STEPS * SECOND_IT_STEPS);
		data_vector data_afm(FIRST_IT_STEPS * SECOND_IT_STEPS);
		data_vector data_sc(FIRST_IT_STEPS * SECOND_IT_STEPS);
		data_vector data_gamma_sc(FIRST_IT_STEPS * SECOND_IT_STEPS);
		data_vector data_xi_sc(FIRST_IT_STEPS * SECOND_IT_STEPS);
		data_vector data_eta(FIRST_IT_STEPS * SECOND_IT_STEPS);

		std::vector<data_vector*> data_mapper = {
			&data_cdw, &data_afm, &data_sc, &data_gamma_sc, &data_xi_sc, &data_eta
		};

		PhaseHelper phaseHelper(input, rank, numberOfRanks);
		phaseHelper.compute_crude(data_mapper);

		data_vector recieve_cdw, recieve_afm, recieve_sc, recieve_gamma_sc, recieve_xi_sc, recieve_eta;
		
		recieve_cdw.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		recieve_afm.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		recieve_sc.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		recieve_gamma_sc.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		recieve_xi_sc.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		recieve_eta.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		
#ifndef _NO_MPI
		MPI_Allgather(data_cdw.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_cdw.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(data_afm.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_afm.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(data_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(data_gamma_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_gamma_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(data_xi_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_xi_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgather(data_eta.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_eta.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, MPI_COMM_WORLD);
#endif
		if (rank == 0) {
			std::vector<std::string> comments;
			comments.push_back(input.getString("second_iterator_type") + "_MIN=" + std::to_string(SECOND_IT_MIN)
				+ "   " + input.getString("second_iterator_type") + "_MAX=" + std::to_string(SECOND_IT_MAX));

			comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
				+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

			std::string output_folder = input.getString("output_folder");
			std::filesystem::create_directories("../../data/phases/" + output_folder);

			Utility::saveData_boost(recieve_cdw, SECOND_IT_STEPS, "../../data/phases/" + output_folder + "cdw.dat.gz", comments);
			Utility::saveData_boost(recieve_afm, SECOND_IT_STEPS, "../../data/phases/" + output_folder + "afm.dat.gz", comments);
			Utility::saveData_boost(recieve_sc, SECOND_IT_STEPS, "../../data/phases/" + output_folder + "sc.dat.gz", comments);
			Utility::saveData_boost(recieve_gamma_sc, SECOND_IT_STEPS, "../../data/phases/" + output_folder + "gamma_sc.dat.gz", comments);
			Utility::saveData_boost(recieve_xi_sc, SECOND_IT_STEPS, "../../data/phases/" + output_folder + "xi_sc.dat.gz", comments);
			Utility::saveData_boost(recieve_eta, SECOND_IT_STEPS, "../../data/phases/" + output_folder + "eta.dat.gz", comments);
		}

		if (input.getBool("improved_boundaries")) {
			data_mapper = {
				&recieve_cdw, &recieve_afm, &recieve_sc, &recieve_gamma_sc, &recieve_xi_sc, &recieve_eta
			};
			int NUMBER_OF_GAP_VALUES = data_mapper.size();

			data_vector local[NUMBER_OF_GAP_VALUES];
			int sizes[NUMBER_OF_GAP_VALUES];
			for (size_t i = 0; i < NUMBER_OF_GAP_VALUES; i++)
			{
				phaseHelper.findSingleBoundary(*(data_mapper[i]), local[i], i, rank);
				sizes[i] = local[i].size();
			}

			std::vector<int> all_sizes[NUMBER_OF_GAP_VALUES];

			data_vector recieve_boundaries[NUMBER_OF_GAP_VALUES];
			if (rank == 0) {
				for (size_t i = 0; i < NUMBER_OF_GAP_VALUES; i++)
				{
					all_sizes[i].resize(numberOfRanks);
				}
			}
			std::vector<int> totalSizes(NUMBER_OF_GAP_VALUES, 0);
#ifndef _NO_MPI
			for (size_t i = 0; i < NUMBER_OF_GAP_VALUES; i++)
			{
				MPI_Gather(&(sizes[i]), 1, MPI_INT, all_sizes[i].data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
				std::vector<int> displacements(numberOfRanks, 0);

				if(rank == 0){
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
				
				std::string names[NUMBER_OF_GAP_VALUES] = {"cdw", "afm", "sc", "gamma_sc", "xi_sc", "eta"};
				for (size_t i = 0; i < NUMBER_OF_GAP_VALUES; i++)
				{
					const int n = recieve_boundaries[i].size() / 2;
					for (size_t j = 1; j < n; j+=2)
					{
						std::swap(recieve_boundaries[i][j], recieve_boundaries[i][n + j - 1]);
					}
					for (int k = 0; k < n; k++) {
    				    for (int l = k + 1; l < n; l++) {
    				        if (recieve_boundaries[i][k] > recieve_boundaries[i][l]) {
    				            std::swap(recieve_boundaries[i][k], recieve_boundaries[i][l]);
    				            std::swap(recieve_boundaries[i][n + k], recieve_boundaries[i][n + l]);
    				        }
    				    }
    				}

					Utility::saveData_boost(recieve_boundaries[i], n, "../../data/phases/" + output_folder + "boundaries_" + names[i] + ".dat.gz", comments);
				}
			}
		}
	}
	else if (input.getString("compute_what") == "modes") {
		omp_set_num_threads(8);

		Hubbard::Model::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
			0, 0, input.getString("global_iterator_type"), input.getString("second_iterator_type"));

		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		std::vector<data_vector> reciever;
		std::vector<data_vector> oneParticleEnergies;
		double totalGapValue;
		std::unique_ptr<std::vector<Hubbard::Resolvent_L>> resolvents;

		std::unique_ptr<Hubbard::Model> model;
		if (input.getBool("use_broyden")) {
			model = std::make_unique<Hubbard::UsingBroyden>(
				Hubbard::UsingBroyden(modelParameters, input.getInt("number_of_basis_terms"), input.getInt("start_basis_at")));
		}
		else {
			model = std::make_unique<Hubbard::HubbardCDW>(
				Hubbard::HubbardCDW(modelParameters, input.getInt("number_of_basis_terms"), input.getInt("start_basis_at")));
		}
		model->computePhases();
		totalGapValue = model->getTotalGapValue();
		resolvents = model->computeCollectiveModes(reciever);
		model->getAllEnergies(oneParticleEnergies);

		if (rank == 0) {
			std::vector<std::string> comments;
			comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
				+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

			std::string output_folder = input.getString("output_folder") + modelParameters.getFileName();
			std::filesystem::create_directories("../../data/" + output_folder);

			comments.push_back("Total Gap=" + std::to_string(totalGapValue));
			if (!(reciever.empty())) {
				Utility::saveData_boost(reciever, "../../data/" + output_folder + ".dat.gz", comments);
			}
			if (resolvents) {
				std::string names[6] = { "phase_SC", "phase_CDW", "phase_AFM", "higgs_SC", "higgs_CDW", "higgs_AFM" };
				for (size_t i = 0; i < resolvents->size(); i++)
				{
					(*resolvents)[i].writeDataToFile("../../data/" + output_folder + "resolvent_" + names[i]);
				}
			}
			else {
				std::cout << "Resolvent returned a null pointer." << std::endl;
			}
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			comments.pop_back();
			Utility::saveData_boost(oneParticleEnergies, "../../data/" + output_folder + "one_particle.dat.gz", comments);
		}
	}

	if(rank == 0){
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	}
#ifndef _NO_MPI
	return MPI_Finalize();
#else
	return 0;
#endif
}