import pysam
import sys,os

variant_type = sys.argv[1]  # "coding" or "noncoding"
chromosome = int(sys.argv[2])

try:
    pysam.tabix_index(f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/{variant_type}/ADSP.chr{chromosome}.vcf", preset="vcf")
except:
    print("Issue with", variant_type, chromosome)