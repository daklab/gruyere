import pandas as pd
import os, sys
import numpy as np
from tqdm import tqdm

PHENOTYPE_FILE = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/phenotypes/36K_QC_filtered_final.csv"
OUTPUT_EXTRACT_PATH = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract/' 
ABC_PATH = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_all/glass_lab_fastq/processed_files/ABC_data"
GENE_POS = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/gene_positions38_ensemble.txt"
OUTPUT_EXTRACT_PATH = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract/' 

def make_extract_coding():
    output = os.path.join(OUTPUT_EXTRACT_PATH, 'coding')
    if not os.path.exists(output): 
        os.mkdir(output)
    gene_pos= pd.read_csv(GENE_POS, sep = "\t")
    gene_pos = gene_pos[gene_pos['biotype']=="protein_coding"]
    for chromosome in range(1, 23):
        outfile = os.path.join(output, f"ADSP.chr{chromosome}.txt")
        chr_pos = gene_pos[gene_pos['chr']==str(chromosome)]
        chr_pos = chr_pos.drop_duplicates("ensgene")
        chr_pos['Gene'] = np.where(chr_pos['symbol'].isna() == False, chr_pos['symbol'], chr_pos['ensgene'])
        chr_pos = chr_pos[['chr','start','end','Gene']].drop_duplicates()
        chr_pos.to_csv(outfile, sep = "\t", header = False, index = False)
    return

def make_extract_noncoding(cell):
    output = os.path.join(OUTPUT_EXTRACT_PATH, cell)
    if not os.path.exists(output): 
        os.mkdir(output)
    cell_df = pd.read_csv(os.path.join(ABC_PATH, f'ABC_results_{cell}/{cell}/Predictions/EnhancerPredictionsFull_threshold0.02_self_promoter.tsv'), sep = "\t")
    cell_df = cell_df[['chr','start','end','TargetGene']]
    cell_df['CHR'] = cell_df['chr'].str.split("chr").str[1]
    for chromosome in range(1, 23):
        outfile = os.path.join(output, f"ADSP.chr{chromosome}.txt")
        chr_df = cell_df[cell_df['CHR']==str(chromosome)]
        chr_df = chr_df[['CHR','start','end','TargetGene']].drop_duplicates()
        chr_df.to_csv(outfile, sep = "\t", header = False, index = False)
    return
        
def make_one_extract_noncoding():
    output = os.path.join(OUTPUT_EXTRACT_PATH, 'noncoding')
    if not os.path.exists(output): 
        os.mkdir(output)
    noncoding_df = pd.DataFrame()
    for cell in ['microglia','oligodendrocyte','astrocyte','neuron']:
        cell_df = pd.read_csv(os.path.join(ABC_PATH, f'ABC_results_{cell}/{cell}/Predictions/EnhancerPredictionsFull_threshold0.02_self_promoter.tsv'), sep = "\t")
        cell_df = cell_df[['chr','start','end','TargetGene']]
        cell_df['CHR'] = cell_df['chr'].str.split("chr").str[1]
        noncoding_df = pd.concat((noncoding_df, cell_df))
        for chromosome in range(1, 23):
            outfile = os.path.join(output, f"ADSP.chr{chromosome}.txt")
            chr_df = noncoding_df[noncoding_df['CHR']==str(chromosome)]
            chr_df = chr_df.drop_duplicates(['CHR','start','end'])
            chr_df = chr_df[['CHR','start','end','TargetGene']].drop_duplicates()
            chr_df.to_csv(outfile, sep = "\t", header = False, index = False)
    return chr_df
    
def main():
    make_extract_coding()
    for cell in ['microglia','oligodendrocyte','astrocyte','neuron']:
        make_extract_noncoding(cell)
    make_one_extract_noncoding()
    return

if __name__ == "__main__":
    main()




