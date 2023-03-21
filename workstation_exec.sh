module purge
module load gcc/11.3.0-full
module load boost/1.71.0
make
mpirun -n 1 ./build/main params/params_modes.config