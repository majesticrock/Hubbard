make -j8
mpirun -n 1 ./build/main params/params_$1.config