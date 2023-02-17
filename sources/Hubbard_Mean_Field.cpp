#define _USE_MATH_DEFINES

#include <omp.h>
#ifdef NDEBUG
#include <mpi.h>
#endif

#include <string>
#include <iostream>
#include <filesystem>
#include <chrono>
#include <memory>

#include "Utility/InputFileReader.hpp"
#include "Utility/OutputWriter.hpp"
#include "Hubbard/BasicHubbardModel.hpp"
#include "Hubbard/HubbardCDW.hpp"
#include "Hubbard/UsingBroyden.hpp"
#include "Hubbard/Constants.hpp"

int Hubbard::Constants::K_DISCRETIZATION = 100;

std::ostream& operator<<(std::ostream& os, const Hubbard::Model::ModelParameters& mp) {
	os << mp.temperature << "\t" << mp.U << "\t" << mp.V;
	return os;
}

int main(int argc, char** argv)
{
#ifdef NDEBUG
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

#ifdef _DO_TEST
	Hubbard::Model::ModelParameters mP(0.1, -2, 0, 0, 0, "", "");
	Hubbard::HubbardCDW model(mP);

	std::chrono::steady_clock::time_point test_b = std::chrono::steady_clock::now();
	model.computePhases(true).print();
	std::chrono::steady_clock::time_point test_e = std::chrono::steady_clock::now();
	std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;

	std::cout << "\n\n";
	model.parseCommutatorData();
	std::cout << "\n\n" << std::endl;
	Hubbard::UsingBroyden model2(mP);

	test_b = std::chrono::steady_clock::now();
	model2.computePhases(true).print();
	test_e = std::chrono::steady_clock::now();
	std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(test_e - test_b).count() << "[ms]" << std::endl;

	std::vector<std::vector<double>> energies;
	model2.getEnergies(energies, 1);
	Utility::saveData(energies, "../../dataenergies.txt");
	return MPI_Finalize();
#endif // _DO_TEST
	// Setup the parameters T, U, V
	std::vector<double> model_params = input.getDoubleList("model_parameters");

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
		data_vector  data_sc(FIRST_IT_STEPS * SECOND_IT_STEPS);
		data_vector data_eta(FIRST_IT_STEPS * SECOND_IT_STEPS);

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
					Hubbard::UsingBroyden model(local);
					ret = model.computePhases();
				}
				else {
					Hubbard::HubbardCDW model(local);
					ret = model.computePhases();
				}

				data_cdw[(T * SECOND_IT_STEPS) + U] = ret.delta_cdw;
				data_sc[(T * SECOND_IT_STEPS) + U] = ret.delta_sc;
				data_eta[(T * SECOND_IT_STEPS) + U] = ret.delta_eta;
				//modelParameters.incrementSecondIterator();
			}
			modelParameters.incrementGlobalIterator();
		}

		std::vector<double> recieve_cdw, recieve_sc, recieve_eta;
		if (rank == 0) {
			recieve_cdw.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
			recieve_sc.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
			recieve_eta.resize(GLOBAL_IT_STEPS * SECOND_IT_STEPS);
		}
#ifdef NDEBUG
		MPI_Gather(data_cdw.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_cdw.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(data_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, recieve_sc.data(), FIRST_IT_STEPS * SECOND_IT_STEPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
			std::filesystem::create_directories("../../data" + output_folder);

			Utility::saveData(recieve_cdw, SECOND_IT_STEPS, "../../data" + output_folder + "cdw.txt", comments);
			Utility::saveData(recieve_sc, SECOND_IT_STEPS, "../../data" + output_folder + "sc.txt", comments);
			Utility::saveData(recieve_eta, SECOND_IT_STEPS, "../../data" + output_folder + "eta.txt", comments);
		}
	}
	else if (input.getString("compute_what") == "modes") {
		omp_set_num_threads(8);
		Hubbard::Model::ModelParameters modelParameters(model_params[0], model_params[1], model_params[2],
			(FIRST_IT_MAX - FIRST_IT_MIN) / FIRST_IT_STEPS, 0,
			input.getString("global_iterator_type"), input.getString("second_iterator_type"));

		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		std::vector<std::vector<data_vector>> reciever(FIRST_IT_STEPS);
		std::vector<std::vector<data_vector>> oneParticleEnergies(FIRST_IT_STEPS);
		std::vector<double> param(FIRST_IT_STEPS);
		std::vector<double> totalGapValues(FIRST_IT_STEPS);

		for (int T = 0; T < FIRST_IT_STEPS; T++)
		{
			std::unique_ptr<Hubbard::Model> model;
			if (input.getBool("use_broyden")) {
				model = std::make_unique<Hubbard::UsingBroyden>(Hubbard::UsingBroyden(modelParameters));
			}
			else {
				model = std::make_unique<Hubbard::HubbardCDW>(Hubbard::HubbardCDW(modelParameters));
			}
			model->computePhases();
			totalGapValues[T] = model->getTotalGapValue();
			model->computeCollectiveModes(reciever[T]);
			model->getAllEnergies(oneParticleEnergies[T]);
			param[T] = modelParameters.getGlobal();
			modelParameters.incrementGlobalIterator();
		}

		if (rank == 0) {
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			std::vector<std::string> comments;
			comments.push_back(input.getString("global_iterator_type") + "_MIN=" + std::to_string(GLOBAL_IT_LIMS[0])
				+ "   " + input.getString("global_iterator_type") + "_MAX=" + std::to_string(GLOBAL_IT_LIMS[1]));

			std::string output_folder = input.getString("output_folder") + input.getString("global_iterator_type") + "_modes/";
			std::filesystem::create_directories("../../data" + output_folder);

			for (int i = 0; i < FIRST_IT_STEPS; i++)
			{
				std::stringstream stream;
				stream << std::fixed << std::setprecision(2) << param[i];
				comments.push_back("Total Gap=" + std::to_string(totalGapValues[i]));
				Utility::saveData(reciever[i], "../../data" + output_folder + stream.str() + ".txt", comments);
				comments.pop_back();
				Utility::saveData(oneParticleEnergies[i], "../../data" + output_folder + stream.str() + "_one_particle.txt", comments);
			}
		}
	}
#ifdef NDEBUG
	return MPI_Finalize();
#else
	return 0;
#endif
}