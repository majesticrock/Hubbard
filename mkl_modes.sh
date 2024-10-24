make MKL -j16
RED='\033[0;31m'
NC='\033[0m' # No Color
n_mpi=1
n_omp=$(($(nproc)/(2*$n_mpi)))
MKL_DYNAMIC=FALSE MKL_NUM_THREADS=8 OMP_NUM_THREADS=$n_omp OMP_MAX_ACTIVE_LEVELS=2 mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build_MKL/HubbardMeanField params/modes.config 2> >(while read line; do echo -e "${RED}$line${NC}"; done) 