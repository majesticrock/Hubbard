make -j8
mpirun -n 4 --map-by node:PE=2 --bind-to core ./build/main params/params_T.config
mpirun -n 4 --map-by node:PE=2 --bind-to core ./build/main params/params_T_finite.config
mpirun -n 4 --map-by node:PE=2 --bind-to core ./build/main params/params_U_positive.config
mpirun -n 4 --map-by node:PE=2 --bind-to core ./build/main params/params_V_positive.config
mpirun -n 4 --map-by node:PE=2 --bind-to core ./build/main params/params_U_negative.config
mpirun -n 4 --map-by node:PE=2 --bind-to core ./build/main params/params_V_negative.config