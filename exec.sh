make -j16
mpirun -n 1 --map-by node:PE=4 --bind-to core ./build/main params/params_$1.config