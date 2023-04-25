make
RED='\033[0;31m'
NC='\033[0m' # No Color
mpirun -n 1 ./build/main params/params_modes.config 2> >(while read line; do echo -e "${RED}$line${NC}"; done) 