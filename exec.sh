make -j8
mpirun -n 4 ./build/main params/params_$1.config