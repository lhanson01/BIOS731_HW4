#!/bin/bash
#SBATCH --array=1-30
#SBATCH --job-name=liam_731_sim
#SBATCH --partition=wrobel
#SBATCH --output=HW4_sim.out
#SBATCH --error=HW4_sim.err

module purge
module load R

JOBID=$SLURM_ARRAY_TASK_ID
Rscript sim/run_simulations.R $JOBID
