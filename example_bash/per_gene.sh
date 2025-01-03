#!/bin/sh

#SBATCH -N 1 
#SBATCH -p pe2
#SBATCH --mem=15G
#SBATCH -t 0-10:00 # Runtime in D-HH:MM
#SBATCH -J gruyere_pergene # <-- name of job
#SBATCH --array=0-21 # <-- number of jobs to run (one per chromosome)

#load required modules
conda activate gruyere

params=("example_data/inputs.yaml")
cell_index=$((SLURM_ARRAY_TASK_ID % ${#params[@]}))
params=${params[$cell_index]}
chr=$((SLURM_ARRAY_TASK_ID / ${#params[@]} + 1))
echo $chr
echo $params

python src/gruyere_pergene.py $params $chr
