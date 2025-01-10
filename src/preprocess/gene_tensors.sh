#!/bin/bash

#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p pe2,bigmem,dev
#SBATCH --mem=200G
#SBATCH -t 2-20:00 # Runtime in D-HH:MM
#SBATCH -J tensors # <-- name of job
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adas@nygenome.org        # Where to send mail
#SBATCH --array=1-22 # <-- number of jobs to run
#SBATCH --output=bash_outputs/stdout_%j.log # Standard output and error log
#SBATCH --error=bash_outputs/error_%j.log

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

module load cuda/10.0
source /gpfs/commons/home/adas/miniconda3/bin/activate
conda activate pyro

python gene_tensors.py $variant_type $chr