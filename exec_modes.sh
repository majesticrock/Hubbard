make -j8
RED='\033[0;31m'
NC='\033[0m' # No Color
mpirun -n 1 --map-by node:PE=8 --bind-to core ./build/main params/params_modes.config 2> >(while read line; do echo -e "${RED}$line${NC}"; done) 