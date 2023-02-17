module purge
module add gcc/11.3.0-full
module add mpi/openmpi3-x86_64
make
mpirun -n 1 ./build/main params/params_modes.config