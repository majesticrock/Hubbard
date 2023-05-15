#!/bin/bash

sbatch slurm/T_finite.slurm
sbatch slurm/T.slurm
sbatch slurm/U_neg.slurm
sbatch slurm/U_pos.slurm
sbatch slurm/V_neg.slurm
sbatch slurm/V_pos.slurm