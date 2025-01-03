#!/bin/bash

#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p pe2,bigmem,dev
#SBATCH --mem=100G
#SBATCH -t 0-20:00 # Runtime in D-HH:MM
#SBATCH -J genotype_plink # <-- name of job
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adas@nygenome.org        # Where to send mail
#SBATCH --array=1-44 # <-- number of jobs to run
#SBATCH --output=bash_outputs/stdout_%j.log # Standard output and error log
#SBATCH --error=bash_outputs/error_%j.log

cells=("coding" "noncoding")
chr=$(( (SLURM_ARRAY_TASK_ID + 1) / 2 ))          # Ceiling of task_id / 2
cell_idx=$(( (SLURM_ARRAY_TASK_ID - 1) % 2 ))     # 0 for "coding", 1 for "noncoding"
cell=${cells[$cell_idx]}              # Get the cell type from the array
echo $chr
echo $cell

keep="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract2/keep_qc.txt"
extract="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract/${cell}/ADSP.chr${chr}.txt"
vcf="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated/ADSP.chr${chr}.vcf.gz"
output="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/${cell}/ADSP.chr${chr}"
echo $vcf
echo $keep
echo $extract
echo $output

module load cuda/10.0
source /gpfs/commons/home/adas/miniconda3/bin/activate
conda activate rv_ad

plink --vcf $vcf --keep $keep  --extract range $extract --out $output --const-fid --maf 0.000000000000001 --recode vcf --biallelic-only --geno 0.1 #--mind 0.1
