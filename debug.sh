make debug -j
n_mpi=1
n_omp=$(($(nproc)/(2*$n_mpi)))
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core gdb ./build/main
# Call this script, then gdb opens.
# Within gdp "run <configfile>", e.g. "run params/params_modes.config"