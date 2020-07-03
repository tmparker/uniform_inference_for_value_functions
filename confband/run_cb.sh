#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=250M
#SBATCH --time=0-00:45 # time (DD-HH:MM)
#SBATCH --array=7,8

module load gcc r
Rscript sim_cband.R $SLURM_ARRAY_TASK_ID

