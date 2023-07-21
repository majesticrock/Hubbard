make -j8
mpirun -n 1 --map-by node:PE=8 --bind-to core ./build/main params/params_test.config