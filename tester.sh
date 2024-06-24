make -j16
n_mpi=1
n_omp=$(($(nproc)/(2*$n_mpi)))
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/HubbardMeanField params/test.config