#!/bin/bash

#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p pe2
#SBATCH --mem=200G
#SBATCH -t 0-20:00 # Runtime in D-HH:MM
#SBATCH -J map_annotations # <-- name of job
#SBATCH --array=1-44 # <-- number of jobs to run
#SBATCH --output=bash_outputs/tstdout_%j.log # Standard output and error log
#SBATCH --error=bash_outputs/terror_%j.log

module load cuda/10.0
source /gpfs/commons/home/adas/miniconda3/bin/activate
conda activate rv_ad

task_id=$SLURM_ARRAY_TASK_ID
if (( task_id <= 22 )); then
    chr=$task_id
    variant_type="coding"
else
    chr=$(( task_id - 22 ))
    variant_type="noncoding"
fi

echo "Chromosome: $chr"
echo "Variant Type: $variant_type"

python map_annotations.py $variant_type $chr


