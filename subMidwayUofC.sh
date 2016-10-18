#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=simu
#SBATCH --mem=32000
#SBATCH --output=to-%A-%a
#SBATCH --array=1
#SBATCH --time=00:04:00

module load R/3.1+intel-15.0

export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE

Rscript proposed.R $SLURM_ARRAY_TASK_ID
