#!/bin/bash

# Define the token to search for and the new value to set
TOKEN="model_parameters"
NEW_VALUES=("0 -2 0.1" "0 -2 0.2" "0 -2 0.3")

for NEW_VALUE in "${NEW_VALUES[@]}"; do
  # Loop through each line in the config file
  while read line; do
    # Skip lines beginning with #
    if [[ $line == \#* ]]; then
      continue
    fi
    
    # Split the line into token and value
    TOKEN_NAME=$(echo "$line" | awk '{print $1}')
    TOKEN_VALUE=$(echo "$line" | cut -d' ' -f2-)
    
    # Check if the current token matches the one we're looking for
    if [[ "$TOKEN_NAME" == "$TOKEN" ]]; then
      # If so, replace the value with the new one
      sed -i "s/$TOKEN_NAME $TOKEN_VALUE/$TOKEN_NAME $NEW_VALUE/" params/params_auto.config
      break
    fi
  done < params/params_auto.config
  
  # Execute the program
  mpirun -n 1 ./build/main params/params_auto.config
  wait
done