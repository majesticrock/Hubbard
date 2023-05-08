#!/bin/bash

make

rm -rf auto_generated/
mkdir -p auto_generated
# Define the name of the input file
input_file="params/for_auto.txt"

# Read the lines of the input file into an array of strings
readarray -t NEW_VALUES < "${input_file}"
TOKEN="model_parameters"


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
      sed "s/$TOKEN_NAME $TOKEN_VALUE/$TOKEN_NAME $NEW_VALUE/" params/params_auto.config > auto_generated/params_$NEW_VALUE.config
      break
    fi
  done < params/params_auto.config
  cp slurm/modes.slurm auto_generated/i.slurm
  sed -i "s/mpirun -n 1 ./build/main params/params_modes.config/mpirun -n 1 ./build/main auto_generated/params_$NEW_VALUE.config" auto_generated/$NEW_VALUE.slurm

  # Execute the program
  sbatch auto_generated/$NEW_VALUE.slurm
done