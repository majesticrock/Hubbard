make -j
RED='\033[0;31m'
NC='\033[0m' # No Color
n_mpi=8
n_omp=$(($(nproc)/(2*$n_mpi)))
mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/main params/unknown_boundary.config 2> >(while read line; do echo -e "${RED}$line${NC}"; done) 