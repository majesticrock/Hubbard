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

	typedef std::vector<double> data_vector;
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

		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		for (int T = 0; T < FIRST_IT_STEPS; T++)
		{
#pragma omp parallel for num_threads(4) schedule(dynamic)
			for (int U = 0; U < SECOND_IT_STEPS; U++)
			{
				Hubbard::Model::ModelParameters local = modelParameters;
				local.setSecondIterator(U);
				Hubbard::Model::data_set ret;
				if (input.getBool("use_broyden")) {
					Hubbard::UsingBroyden model(local, 0, 0);
					ret = model.computePhases();
				}
				else {
					Hubbard::HubbardCDW model(local, 0, 0);
					ret = model.computePhases();
				}

				data_cdw[(T * SECOND_IT_STEPS) + U] = ret.delta_cdw;
				data_afm[(T * SECOND_IT_STEPS) + U] = ret.delta_afm;
				data_sc[(T * SECOND_IT_STEPS) + U] = ret.delta_sc;
				data_gamma_sc[(T * SECOND_IT_STEPS) + U] = ret.gamma_sc;
				data_xi_sc[(T * SECOND_IT_STEPS) + U] = ret.xi_sc;
				data_eta[(T * SECOND_IT_STEPS) + U] = ret.delta_eta;
				//modelParameters.incrementSecondIterator();
			}
			modelParameters.incrementGlobalIterator();
		}
		if (input.getBool("improved_boundaries")) {
			constexpr int NUMBER_OF_GAP_VALUES = 6;

			modelParameters = Hubbard::Model::ModelParameters(model_params[0], model_params[1], model_params[2],
				(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, 0,
				input.getString("global_iterator_type"), input.getString("second_iterator_type"));

			data_vector local_data(2 * NUMBER_OF_GAP_VALUES * FIRST_IT_STEPS);

			for (int T = 0; T < FIRST_IT_STEPS; T++)
			{
				for (size_t i = 0; i < NUMBER_OF_GAP_VALUES; i++)
				{
					double lower_boundary = -128;
					double upper_boundary = 128;
					bool found_lower = false, found_upper = false;

					for (size_t j = 1; j < SECOND_IT_STEPS - 1; j++)
					{
						if (std::abs( (*(data_mapper[i]))[T * SECOND_IT_STEPS + j] ) > 1e-12) {
							if (std::abs( (*(data_mapper[i]))[T * SECOND_IT_STEPS + j - 1] ) < 1e-12) {
								// Current element is nonzero and the previous element is 0
								lower_boundary = modelParameters.setSecondIterator(j);
								if (i == 0) {
									std::cout << "cdw:  " << lower_boundary << std::endl;
								}
								found_lower = true;
							}
							if (std::abs( (*(data_mapper[i]))[T * SECOND_IT_STEPS + j + 1] ) < 1e-12) {
								// Current element is nonzero and the next element is 0
								upper_boundary = modelParameters.setSecondIterator(j);
								found_upper = true;
							}
						}
					}
					if (!found_lower) {
						if (std::abs((*data_mapper[i])[T * SECOND_IT_STEPS]) > 1e-12) {
							// Lowest element is nonzero, thus the lower boundary, if it exists,
							// is outside the scope of the current simulation
							local_data[2 * i + FIRST_IT_STEPS * T] = -std::numeric_limits<double>::infinity();
						}
						else {
							// The lowest element is zero, thereby all higher elements are 0 too
							// This parameter is never nonzero.
							local_data[2 * i + FIRST_IT_STEPS * T] = std::numeric_limits<double>::quiet_NaN();
							local_data[1 + 2 * i + FIRST_IT_STEPS * T] = std::numeric_limits<double>::quiet_NaN();
							break;
						}
					}
					if (!found_upper) {
						if (std::abs((*data_mapper[i])[(T + 1) * SECOND_IT_STEPS - 1]) > 1e-12) {
							// Highest element is nonzero, thus the upper boundary, if it exists,
							// is outside the scope of the current simulation
							local_data[1 + 2 * i + FIRST_IT_STEPS * T] = std::numeric_limits<double>::infinity();
						}
						// The other case would imply that all elements are 0
						// This is already covered in the previous if-statement
					}

					double convergence;
					double current_lower;
					double current_upper;
					double current_test;
					if (found_lower) {
						std::cout << lower_boundary << std::endl;
						current_lower = lower_boundary - modelParameters.getSecondStep();
						current_upper = lower_boundary;
						while ((convergence = (current_upper - current_lower)) > 1e-9) {
							current_test = current_lower + 0.5 * convergence;

							modelParameters.setSecondIteratorExact(current_test);
							Hubbard::Model::data_set ret;
							if (input.getBool("use_broyden")) {
								Hubbard::UsingBroyden model(modelParameters, 0, 0);
								ret = model.computePhases();
							}
							else {
								Hubbard::HubbardCDW model(modelParameters, 0, 0);
								ret = model.computePhases();
							}

							if (ret.isFinite(i)) {
								current_upper = current_test;
							}
							else {
								current_lower = current_test;
							}
						}
						local_data[2 * i + FIRST_IT_STEPS * T] = current_upper;
					}

					if (found_upper) {
						std::cout << upper_boundary << std::endl;
						current_lower = upper_boundary;
						current_upper = upper_boundary + modelParameters.getSecondStep();
						while ((convergence = (current_upper - current_lower)) > 1e-8) {
							current_test = current_lower + 0.5 * convergence;

							modelParameters.setSecondIteratorExact(current_test);
							Hubbard::Model::data_set ret;
							if (input.getBool("use_broyden")) {
								Hubbard::UsingBroyden model(modelParameters, 0, 0);
								ret = model.computePhases();
							}
							else {
								Hubbard::HubbardCDW model(modelParameters, 0, 0);
								ret = model.computePhases();
							}

							if (ret.isFinite(i)) {
								current_lower = current_test;
							}
							else {
								current_upper = current_test;
							}
						}
						local_data[1 + 2 * i + FIRST_IT_STEPS * T] = current_lower;
					}
				}
				modelParameters.incrementGlobalIterator();
			}

			data_vector recieve_boundaries;
			if (rank == 0) {
				recieve_boundaries.resize(2 * NUMBER_OF_GAP_VALUES * GLOBAL_IT_STEPS);
			}
#ifndef _NO_MPI
			MPI_Gather(local_data.data(), 2 * NUMBER_OF_GAP_VALUES * FIRST_IT_STEPS, MPI_DOUBLE, recieve_boundaries.data(), 2 * NUMBER_OF_GAP_VALUES * FIRST_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
			if (rank == 0) {
				std::vector<std::string> comments;
				comments.push_back(input.getString("second_iterator_type") + "_MIN=" + std::to_string(SECOND_IT_MIN)
					+ "   " + input.getString("second_iterator_type") + "_MAX=" + std::to_string(SECOND_IT_MAX));

				comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
					+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

				std::string output_folder = input.getString("output_folder");
				std::filesystem::create_directories("../../data/phases/" + output_folder);

				Utility::saveData_boost(recieve_boundaries, SECOND_IT_STEPS, "../../data/phases/" + output_folder + "boundaries.dat.gz", comments);
			}
		}

		data_vector recieve_cdw, recieve_afm, recieve_sc, recieve_gamma_sc, recieve_xi_sc, recieve_eta;
		if (rank == 0) {
			recieve_cdw.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
			recieve_afm.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
			recieve_sc.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
			recieve_gamma_sc.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
			recieve_xi_sc.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
			recieve_eta.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		}
#ifndef _NO_MPI
		MPI_Gather(data_cdw.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_cdw.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(data_afm.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_afm.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(data_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(data_gamma_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_gamma_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(data_xi_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_xi_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(data_eta.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_eta.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
		if (rank == 0) {
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

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
	}
	else if (input.getString("compute_what") == "modes") {
		omp_set_num_threads(8);

		Hubbard::Model::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
			(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, 0,
			input.getString("global_iterator_type"), input.getString("second_iterator_type"));

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
#ifndef _NO_MPI
	return MPI_Finalize();
#else
	return 0;
#endif
}