#ifdef _DEBUG
#define _NO_MPI
#endif
#include "Hubbard/GlobalDefinitions.hpp"

#include <omp.h>
#ifndef _NO_MPI
#define _DEFAULT_EXIT MPI_Finalize()
#include <mpi.h>
#else
#define _DEFAULT_EXIT 0
#endif

#include <chrono>
#include "Hubbard/Constants.hpp"
#include "TestHandler.hpp"
#include "PhaseHandler.hpp"
#include "ModeHandler.hpp"

using namespace Hubbard;

int Constants::K_DISCRETIZATION = -1;
int Constants::BASIS_SIZE = -1;
int Constants::HALF_BASIS = -1;
global_floating_type Constants::PI_DIV_DISCRETIZATION = -1;

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
	Constants::setDiscretization(input.getInt("k_discretization"));

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