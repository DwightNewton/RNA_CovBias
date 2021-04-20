#!/bin/sh
#SBATCH --time=06:00:00
#SBATCH --array=1-379
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END

sed "${SLURM_ARRAY_TASK_ID}q;d" picard_mdd.cmdlist | bash