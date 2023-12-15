#!/bin/bash

## For more information on the following SBATCH commands, please refer to the manual page (man sbatch) in the terminal
## or look up the slurm documentation under https://slurm.schedmed.com/sbatch.html

## Mandatory:
#SBATCH --job-name=small_modes
#SBATCH --output=/home/althueser/phd/cpp/HubbardMeanField/small_modes_output_0_4.275_1.txt
#SBATCH --time=16:00:00         ## maximum runtime; hours:minutes:seconds
#SBATCH --partition=ultralong       ## choose queue

#SBATCH --ntasks=1             ## sets number of total mpi processes
#SBATCH --tasks-per-node=1    ## mpi processes spawned per node
#SBATCH --cpus-per-task=8       ## logical cores used per mpi process: should be 1, except if you want to combine OpenMP with$

#SBATCH --mem=36gb               ## sets maximum allowed memory per node
##SBATCH --mem-per-cpu=100mb    ## sets maximum allowed memory per mpi process

#SBATCH --mail-user=joshua.althueser@tu-dortmund.de     ## replace mail by personal mail address
#SBATCH --mail-type=FAIL ## most relevant options: NONE, BEGIN, END, FAIL

## Optional:
#SBATCH --hint=nomultithread    ## deactivate Hyperthreading (recommended); for Hyperthreading comment out this line
#SBATCH --constraint=CascadeLake       ## chose a specific feature, e.g., only nodes with Haswell-architecture
                                ## Feature-Output by "cat /etc/slurm/slurm.conf | grep Feature"
module purge
module add gcc/11.3.0-full
module add boost/1.71.0
cd /home/althueser/phd/cpp/HubbardMeanField/        ## go to working directory

date
echo "--- START ---"

## execute binary using mpirun on the allocated computation resource; the number of cores is $
mpirun ./build/main small_auto_generated2/0_4.275_1.config

echo "--- END ---"
date
echo
echo "$(whoami) is leaving from $(hostname) ..."
echo
