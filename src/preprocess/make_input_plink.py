import pandas as pd
import os, sys
import numpy as np
from tqdm import tqdm


variant_type = sys.argv[1]  # "coding" or "noncoding"
chromosome = int(sys.argv[2])

PHENOTYPE_FILE = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/phenotypes/36K_QC_filtered_final.csv"
OUTPUT_EXTRACT_PATH = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract/' 
ABC_PATH = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_all/glass_lab_fastq/processed_files/ABC_data"
GENE_POS = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/gene_positions38_ensemble.txt"
OUTPUT_EXTRACT_PATH = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract/' 
WGSA_PATH = f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/wgsa_annotated/"
GWAS_SIGNIFICANT = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/LD_bellenguez/bellenguez_significant.tsv"

def make_extract(gwas):
    output = os.path.join(OUTPUT_EXTRACT_PATH, variant_type)
    wgsa = pd.read_csv(os.path.join(WGSA_PATH, variant_type, f"WGSA_chr{chromosome}.csv"), usecols = ['chr','pos'])
    wgsa['label'] = "X"
    wgsa['chr'] = wgsa['chr'].astype(int)
    wgsa['pos'] = wgsa['pos'].astype(int)
    wgsa = wgsa[['chr', 'pos','pos', 'label']]
    wgsa = pd.concat((wgsa, gwas))
    outfile = os.path.join(output, f"ADSP.chr{chromosome}.txt")
    wgsa.to_csv(outfile, sep = "\t", header = False, index = False)
    return

def load_GWAS_sig():
    output = os.path.join(OUTPUT_EXTRACT_PATH, variant_type)
    if not os.path.exists(output): 
        os.mkdir(output)
    gwas = pd.read_csv(GWAS_SIGNIFICANT)
    gwas = gwas[gwas['CHR']==f"chr{chromosome}"]
    outfile = os.path.join(output, f"Bellenguez_significant.chr{chromosome}.txt")
    gwas[['SNP']].to_csv(outfile, header = False, index = False)
    gwas['chr'] = chromosome
    gwas['pos'] = gwas['BP'].astype(int)
    gwas['label'] = "GWAS"
    gwas = gwas[['chr', 'pos','pos', 'label']]
    return gwas

    
def main():
    gwas = load_GWAS_sig()
    make_extract(gwas)
    return
    
if __name__ == "__main__":
    main()