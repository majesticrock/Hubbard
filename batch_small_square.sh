#!/bin/bash
rm -rf small_square_auto_generated/
mkdir -p small_square_auto_generated
# Define the name of the input file
input_file="params/for_auto_square.txt"

# Read the lines of the input file into an array of strings
readarray -t NEW_VALUES < "${input_file}"
TOKEN="model_parameters"


for NEW_VALUE in "${NEW_VALUES[@]}"; do
  NEW_NAME=$(echo "$NEW_VALUE" | sed 's/ /_/g')
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
      sed "s/$TOKEN_NAME $TOKEN_VALUE/$TOKEN_NAME $NEW_VALUE/" params/small_square.config > small_square_auto_generated/$NEW_NAME.config
      break
    fi
  done < params/small_square.config
  cp slurm/modes_small_square.slurm small_square_auto_generated/$NEW_NAME.slurm
  sed -i "s|#SBATCH --job-name=modes|#SBATCH --job-name=$NEW_NAME|" small_square_auto_generated/$NEW_NAME.slurm
  sed -i "s|#SBATCH --output=/home/althueser/phd/cpp/HubbardMeanField/small_modes_output.txt|#SBATCH --output=/home/althueser/phd/cpp/HubbardMeanField/small_square_output_$NEW_NAME.txt|" small_square_auto_generated/$NEW_NAME.slurm
  sed -i "s|mpirun ./build/HubbardMeanField params/small_square.config|mpirun ./build/HubbardMeanField small_square_auto_generated/$NEW_NAME.config|" small_square_auto_generated/$NEW_NAME.slurm

  # Execute the program
  sbatch small_square_auto_generated/$NEW_NAME.slurm
done