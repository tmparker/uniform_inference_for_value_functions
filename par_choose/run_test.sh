#!/bin/sh
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=5M
#SBATCH --time=00:01
#SBATCH --array=1-2

module load r
echo "testing $SLURM_ARRAY_TASK_ID"
Rscript test.R $SLURM_ARRAY_TASK_ID


