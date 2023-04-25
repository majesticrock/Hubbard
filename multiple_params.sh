#!/bin/bash

make
# Define the token to search for and the new value to set
TOKEN="model_parameters"
NEW_VALUES=("0 -2 -2" "0 -2 -1.9" "0 -2 -1" "0 -2 -0.1")
#NEW_VALUES=("0 -2 -0.1")

RED='\033[0;31m'
NC='\033[0m' # No Color

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
      sed -i "s/$TOKEN_NAME $TOKEN_VALUE/$TOKEN_NAME $NEW_VALUE/" params/params_auto.config
      break
    fi
  done < params/params_auto.config
  
  # Execute the program
  mpirun -n 1 ./build/main params/params_auto.config 2> >(while read line; do echo -e "${RED}$line${NC}"; done) 
  wait
done