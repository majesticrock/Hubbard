#!/bin/bash

sbatch slurm/params_T_finite.config
sbatch slurm/params_T.config
sbatch slurm/params_U_negative.config
sbatch slurm/params_U_positive.config
sbatch slurm/params_V_negative.config
sbatch slurm/params_V_positive.config