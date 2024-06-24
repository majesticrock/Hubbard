#!/bin/bash

make -j

# Define the name of the input file
input_file="params/for_auto.txt"

# Read the lines of the input file into an array of strings
readarray -t NEW_VALUES < "${input_file}"
TOKEN="model_parameters"

RED='\033[0;31m'
NC='\033[0m' # No Color

n_mpi=1
n_omp=$(($(nproc)/(2*$n_mpi)))

for NEW_VALUE in "${NEW_VALUES[@]}"; do
  # Loop through each line in the config file
  while read line; do
    if [[ $line == \#* ]]; then
      continue
    fi
    
    # Split the line into token and value
    TOKEN_NAME=$(echo "$line" | awk '{print $1}')
    TOKEN_VALUE=$(echo "$line" | cut -d' ' -f2-)
    
    if [[ "$TOKEN_NAME" == "$TOKEN" ]]; then
      # replace the value with the new one
      sed -i "s/$TOKEN_NAME $TOKEN_VALUE/$TOKEN_NAME $NEW_VALUE/" params/auto.config
      break
    fi
  done < params/auto.config
  
  # Execute the program
  mpirun -n $n_mpi --map-by node:PE=$n_omp --bind-to core ./build/HubbardMeanField params/auto.config 2> >(while read line; do echo -e "${RED}$line${NC}"; done)
  wait
done