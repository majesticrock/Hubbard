make -j8
n_mpi=1
n_omp=$(($(nproc)/(2*$n_mpi)))
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/params_T.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/params_T_finite.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/params_U_positive.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/params_V_positive.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/params_U_negative.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/params_V_negative.config