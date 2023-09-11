make -j16
mpirun -n 1 gdb ./build/main
# Call this script, then gdb opens.
# Within gdp "run <configfile>", e.g. "run params/params_modes.config"