make -j
n_mpi=1
n_omp=$(($(nproc)/(2*$n_mpi)))
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/T.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/T_finite.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/U_positive.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/V_positive.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/U_negative.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/V_negative.config