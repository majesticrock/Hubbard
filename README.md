# Hubbard

The code allows studying the phase diagram and collective excitations of the extended Hubbard model on a square and a simple cubic lattice.
The code was used to compute the mean-field data and collective excitation spectra presented in 

Collective excitations in competing phases in two and three dimensions
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.109.205153?ft=1
doi: https://doi.org/10.1103/PhysRevB.109.205153



## Parameters

The parameter files in the `params` directory are used to control the system parameters.

### select which computation should be performed 
`compute_what <phases|modes|unknown_boundary|test>`
### Use the improved, but more demanding algorithm for the phase diagram
`improved_boundaries <bool>`
### Setup the extended Hubbard model with T, U, V
`model_parameters <float> <float> <float>`
### Allows favouring CDW over SC, default 0.5, 0 = pure SC, 1 = pure CDW
`ratio_CDW_SC <float>`
### Which lattice?
`lattice_type <chain|square|cube>`
### Abandoned project part, would have tried to include electromagnetic coupling
`em_coupling <bool>`
### Whether to use the DOS formulation (see the article) or plain momentum sums
`use_DOS <bool>`
### Number of discretization points for the iterator that is being split to MPI (used for the phase diagram)
`global_iterator_steps <int>`
### Upper end for the MPI iteration
`global_iterator_upper_limit <float>`
### Which parameter should be changed by MPI
`global_iterator_type <U|V|T>`
### Number of discretization points for the second iterator (used for the phase diagram)
`second_iterator_steps <int>`
### Upper end for the inner iteration
`second_iterator_upper_limit <float>`
### Which parameter should be changed by the inner iteration
`second_iterator_type <U|V|T>`
### Half the number of points in k-space per direction
`k_discretization <int>`
### Data dir
`output_folder <string>/`
### Whether to use Broyden's method to solve the self-consistency
`use_broyden <bool>`
### Allows restricting the iEoM basis. If -1, then the full basis is used
`start_basis_at <int>`
### Number of basis terms - default: 12
`number_of_basis_terms <int>`



### Required externals
- Eigen https://eigen.tuxfamily.org/index.php?title=Main_Page or https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html
- Boost https://www.boost.org/
- OpenMPI https://www.open-mpi.org/
- OpenMP https://www.openmp.org/
- CMake https://cmake.org/
- nlohmann/json.hpp https://github.com/nlohmann/json

It should be easy to include MKL, though that has not happened yet (compare https://github.com/majesticrock/ContinuumSystem)
- Intel MKL https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html 

### Required internals

- mrock/utility
- mrock/symbolic_operators

see https://github.com/majesticrock/PhdUtility


## Building

Build with `make`.
The Makefile handles the calls to cmake for you.
Without specification the project will be built for the local machine.
Specifying `cluster` will built for the CPU architecture on the compute cluster.
Could be extended to specify IceLake or CascadeLake, see https://github.com/majesticrock/ContinuumSystem 

## Running the program

Before executing the program, make sure to build and run FermionCommute
https://github.com/majesticrock/FermionCommute/
such that you have the directory `../commutators/hubbard/` filled with the results of the commutators.
`./exec.sh` will execute the program with the parameter file `params/params.config`.
For large scale computations, SLURM scripts are provided.
A few additional bash scripts are provided for ease of running multiple jobs.
