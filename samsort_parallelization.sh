#!/bin/sh
#SBATCH --time=06:00:00
#SBATCH --array=1-379
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END

sed "${SLURM_ARRAY_TASK_ID}q;d" sort_script.cmdlist | bash
