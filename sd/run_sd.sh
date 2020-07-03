#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=250M
#SBATCH --time=00-01:15 # time (DD-HH:MM)
#SBATCH --array=34,45-48,51,52,55,57,77,78,81-88

module load gcc r
Rscript sim_sd.R $SLURM_ARRAY_TASK_ID

