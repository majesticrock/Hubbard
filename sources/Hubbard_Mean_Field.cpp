#ifdef _DEBUG
#define _NO_MPI
#endif

#include <omp.h>
#ifndef _NO_MPI
#define _DEFAULT_EXIT MPI_Finalize()
#include <mpi.h>
#else
#define _DEFAULT_EXIT 0
#endif

#include <chrono>
#include <iostream>
#include "Hubbard/Constants.hpp"
#include "TestHandler.hpp"
#include "PhaseHandler.hpp"
#include "ModeHandler.hpp"

int Hubbard::Constants::K_DISCRETIZATION = 100;
int Hubbard::Constants::BASIS_SIZE = 10000;
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

	if (input.getString("compute_what") == "test") {
		TestHandler test(input, rank, numberOfRanks);
		test.execute(input);
	}
	else if (input.getString("compute_what") == "phases") {
		PhaseHandler phases(input, rank, numberOfRanks);
		phases.execute(input);
	}
	else if (input.getString("compute_what") == "modes") {
		ModeHandler modes(input, rank, numberOfRanks);
		modes.execute(input);
	}

	if (rank == 0) {
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	}
	return _DEFAULT_EXIT;
}