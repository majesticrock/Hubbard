make -j8
mpirun -n 4 --map-by node:PE=2 --bind-to core ./build/main params/params_$1.config