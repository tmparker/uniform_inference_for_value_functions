#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=500M
#SBATCH --time=00-02:30 # time (DD-HH:MM)
#SBATCH --array=3

module load gcc r
Rscript sim_tune.R $SLURM_ARRAY_TASK_ID

