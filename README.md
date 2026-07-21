# Hubbard

The code allows studying the phase diagram and collective excitations of the extended Hubbard model on a square and a simple cubic lattice.
The code was used to compute the mean-field data and collective excitation spectra presented in 

Collective excitations in competing phases in two and three dimensions
J. Althüser and G. S. Uhrig
https://doi.org/10.1103/PhysRevB.109.205153



### Requirements
- C++ 20 and a functioning compiler (tested with g++ 13.3.0 on WSL and g++ 11.5.0 on Red Hat)
- Eigen (tested with version 3.4.1) https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html
- Boost (tested with version 1.78.0) https://www.boost.org/
- OpenMPI (tested with version 4.1.1) https://www.open-mpi.org/
- OpenMP https://www.openmp.org/
- nlohmann/json.hpp (tested with version 3.11.3) https://github.com/nlohmann/json
- [recommended] CMake 3.30 or newer (tested with version 3.31.8)  https://cmake.org/
- [optional] Intel MKL for BLAS and LAPACK (tested with version 2025.2.1) https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html 

- the mrock library located in `../../PhdUtility/`.
- Before executing the program, make sure to build and run `FermionCommute` in `../FermionCommute/`.



## Parameters

The parameter files in the `params` directory are used to control the system parameters.

| Option | Type / Values | Description |
|---|---:|---|
| `compute_what` | `phases` \| `modes` \| `unknown_boundary` \| `test` | Select which computation should be performed. |
| `improved_boundaries` | `<bool>` | Use the improved (but more demanding) algorithm for the phase diagram. |
| `model_parameters` | `<float> <float> <float>` | Setup the extended Hubbard model with `T`, `U`, `V` (in that order). |
| `ratio_CDW_SC` | `<float>` | Allows favouring CDW over SC; default `0.5`. `0` = pure SC, `1` = pure CDW. |
| `lattice_type` | `square` \| `cube` | Which lattice. |
| `use_DOS` | `<bool>` | Whether to use the DOS formulation (see the article) or plain momentum sums. |
| `global_iterator_steps` | `<int>` | Number of discretization points for the MPI-split iterator (used for the phase diagram). |
| `global_iterator_upper_limit` | `<float>` | Upper end for the MPI iteration. |
| `global_iterator_type` | `U` \| `V` \| `T` | Which parameter is changed by MPI. |
| `second_iterator_steps` | `<int>` | Number of discretization points for the inner iterator (used for the phase diagram). |
| `second_iterator_upper_limit` | `<float>` | Upper end for the inner iteration. |
| `second_iterator_type` | `U` \| `V` \| `T` | Which parameter is changed by the inner iteration. |
| `discretization` | `<int>` | Half the number of points in k-space per direction (or entire gamma-space if DOS is used). |
| `output_folder` | `<string>/` | Data directory (note trailing slash). |
| `use_broyden` | `<bool>` | Whether to use Broyden's method to solve the self-consistency. |
| `start_basis_at` | `<int>` | Allows restricting the iEoM basis. `-1` = use full basis. |
| `number_of_basis_terms` | `<int>` | Number of basis terms (default: `12`). |


## Building

Build with `make`.
The Makefile handles the calls to cmake for you.
Without specification the project will be built for the local machine (-march=native).
There are additionally the targets `icelake` and `cascadelake`, which will built for the corresponding CPU architecture.

## Running the program

Before executing the program, make sure to build and run `FermionCommute` in `../FermionCommute/`.
such that you have the directory `../commutators/hubbard/` filled with the results of the commutators.
Then, (after building of course), you may create or edit a parameter file in the `params` directory.
Run the program with `./path/to/executable path/to/param/file.config`.
By default, the executable will be located in `./build/default/hubbard`.
The result will be saved into `../../data/hubbard/<square/cube>/<output_folder in the parameter file>/<modes/phase>/<string generated from the model parameters>/<filename (e.g. resolvents.json.gz)>`.

## Testing

`make test` will build and run the program with a set of parameters listed in the `tests` directory.
It will then call the plot script in the same dir to visualize the test data.
As before, it requires a previously completed run of `FermionCommute`, which must save the commutators to `../commutators/hubbard/`.
The program will save the simulation data and consequently plot it.

The plots will be saved to `build/default/<plot_name>.pdf`.
For doing so, python with matplotlib, numpy, and pandas is required.
Moreover, the python modules of `../../PhdUtility/python` are needed to evaluate the continued fractions.

### Expected results
#### SC-CDW
Peak at ω=0 in Phase, Higgs, and CDW
Peak at ω=2Δ in Higgs and CDW
Both AFM are only inside the continuum

#### SC
Peak at ω=0 in Phase
Peak at ω=2Δ in Higgs
Weak peak at ω<2Δ in CDW
Both AFM are only inside the continuum

#### CDW
Peak at ω=2Δ in Higgs and Phase
Peak at ω=2Δ in CDW
Both AFM are only inside the continuum

#### AFM
SC and CDW only in Continuum
Peak at ω=0 in t.AFM
Peak at ω=2Δ in l.AFM 