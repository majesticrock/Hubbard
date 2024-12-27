#ifdef _DEBUG
#define _NO_MPI
#endif

#ifndef _NO_MPI
#define _DEFAULT_EXIT MPI_Finalize()
#include <mpi.h>
#else
#define _DEFAULT_EXIT 0
#endif
#include <chrono>
#include <Eigen/Dense>
#include "Hubbard/GlobalDefinitions.hpp"
#include "Hubbard/Constants.hpp"
#include "Handler/TestHandler.hpp"
#include "Handler/PhaseHandler.hpp"
#include "Handler/ModeHandler.hpp"
#include "Handler/UnknownBoundaryHandler.hpp"
#include "Handler/ModeDispersionHandler.hpp"

int Hubbard::Constants::K_DISCRETIZATION = -1;
int Hubbard::Constants::BASIS_SIZE = -1;
int Hubbard::Constants::HALF_BASIS = -1;
int Hubbard::Constants::QUARTER_BASIS = -1;
int Hubbard::Constants::EIGHTH_BASIS = -1;
Hubbard::global_floating_type Hubbard::Constants::PI_DIV_DISCRETIZATION = -1;

int main(int argc, char** argv)
{
	using namespace Hubbard;
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
	mrock::Utility::InputFileReader input(argv[1]);
	Constants::setDiscretization(input.getInt("k_discretization"));
	if (input.getBool("use_DOS")) {
		Constants::setBasis(Constants::K_DISCRETIZATION);
	}
	else {
		if (input.getString("lattice_type") == "square") {
			Constants::setBasis(4 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
		}
		else if (input.getString("lattice_type") == "cube") {
			Constants::setBasis(8 * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION * Constants::K_DISCRETIZATION);
		}
		else {
			Constants::setBasis(2 * Constants::K_DISCRETIZATION);
		}
	}
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
	else if (input.getString("compute_what") == "unknown_boundary") {
		UnknownBoundaryHandler u_boundary(input, rank, numberOfRanks);
		u_boundary.execute(input);
	}
	else if (input.getString("compute_what") == "dispersions") {
		assert(argc > 2);
		ModeDispersionHandler dispersions(input, rank, numberOfRanks, std::stoi(argv[2]));
		dispersions.execute(input);
	}

	if (rank == 0) {
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Total runtime = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
	}
	return _DEFAULT_EXIT;
}