make -j16
n_mpi=1
n_omp=$(($(nproc)/(2*$n_mpi)))
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/HubbardMeanField params/T.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/HubbardMeanField params/T_finite.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/HubbardMeanField params/U_positive.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/HubbardMeanField params/V_positive.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/HubbardMeanField params/U_negative.config
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/HubbardMeanField params/V_negative.config