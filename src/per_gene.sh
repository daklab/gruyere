#!/bin/sh
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p pe2
#SBATCH --mem=10G
#SBATCH -t 0-7:00 # Runtime in D-HH:MM
#SBATCH -J parmigiano # <-- name of job
#SBATCH --array=1-2199 #2199 # <-- number of jobs to run 
#SBATCH --output=bash_outputs/tstdout_%j.log               # Standard output and error log
#SBATCH --error=bash_outputs/terror_%j.log
#SBATCH --mail-user=adas@nygenome.org 


#load required modules
module load cuda/10.0
source /gpfs/commons/home/adas/miniconda3/bin/activate
conda activate pyro

params="inputs.yaml"
chr_array=({1..22})
chr_index=$(( SLURM_ARRAY_TASK_ID % 22 ))      # Get the chromosome index (cycles every 22 jobs)
iteration=$(( SLURM_ARRAY_TASK_ID / 22 + 1 ))  # Calculate the iteration (increments every 22 jobs)

chr=${chr_array[$chr_index]}    # Select chromosome from array

echo "Chromosome: $chr"
echo "Iteration: $iteration"
echo "Params: $params"

python per_gene.py $params $chr $iteration
