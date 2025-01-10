# This script generates tensors and gene mappings:
import pysam
import torch
import numpy as np
import pandas as pd
import sys
import os
from tqdm import tqdm
import pickle


variant_type = sys.argv[1]  # "coding" or a cell type
chromosome = int(sys.argv[2])

VARIANTS_EXTRACT = f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract/{variant_type}/ADSP.chr{chromosome}.txt"
ABC_PATH = "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_all/glass_lab_fastq/processed_files/ABC_data"
VCF_PATH = f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/{variant_type}/ADSP.chr{chromosome}.vcf"
WGSA_PATH = f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/wgsa_annotated/{variant_type}/WGSA_chr{chromosome}.csv"
OUTPUT_PATH = f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/gene_matrices/{variant_type}/"

def load_ABC():
    cell_df = pd.DataFrame()
    for cell in ['microglia','oligodendrocyte','astrocyte','neuron']:
        temp = pd.read_csv(os.path.join(ABC_PATH, f'ABC_results_{cell}/{cell}/Predictions/EnhancerPredictionsFull_threshold0.02_self_promoter.tsv'), sep = "\t")
        temp['cell'] = cell
        cell_df = pd.concat((cell_df, temp))
    cell_df = cell_df[cell_df['chr'] == 'chr' + str(chromosome)]
    return cell_df


def load_gene(wgsa):
    vcf = pysam.VariantFile(VCF_PATH)
    genotype_array = []
    variants = []
    for pos in tqdm(list(wgsa['pos'])):
        records = vcf.fetch(str(chromosome), pos - 1, pos)  # VCF is 0-based
        for record in records:
            if record.pos == pos:
                gts = [s['GT'] for s in record.samples.values()]
                if gts:  # Ensure gts is not empty
                    variants.append(pos)
                    genotype_array.append(gts)
                    break
    genotype_array = np.array(genotype_array)  # Ensures it handles None values.
    genotype_array = np.where(genotype_array == None, np.nan, genotype_array).astype(float)
    if genotype_array.ndim == 3:  # Ensure it's 3-dimensional
        genotype_array = np.sum(genotype_array, axis=2)  # Sum alleles for each genotype
        wgsa = wgsa[wgsa['pos'].isin(variants)]
        wgsa['MAF'] = np.nansum(genotype_array, axis=1) / (2 * np.sum(~np.isnan(genotype_array), axis=1))
        variants = list(wgsa['variant_id'])
        return np.array(wgsa.set_index(['variant_id','chr'])), genotype_array, variants
    else:
        print("Issue with genotype array")
        return None

def load_genes(extract, wgsa):
    annot = []
    genotypes = []
    variants = []
    genes = []
    for gene_idx in tqdm(range(extract.shape[0])):
        print("gene index", gene_idx)
        wgsa_sub = wgsa[(wgsa['pos']>= extract.loc[gene_idx,'start']) & (wgsa['pos']<= extract.loc[gene_idx,'end'])]
        try:
            w, g, v = load_gene(wgsa_sub)
        except:
            print("Issue with gene:", gene_idx)
            continue
        annot.append(w)
        genotypes.append(g)
        variants.append(v)
        genes.append([extract.loc[gene_idx,'gene']] * len(v))
    return annot, genotypes, variants, genes


def save_outputs(annot, genotypes, variants, genes):
    df = pd.DataFrame(np.column_stack((np.hstack(genes), np.hstack(variants))))
    df = pd.DataFrame(np.column_stack((np.hstack(genes), np.hstack(variants))))
    gene_to_index = {gene: idx for idx, gene in enumerate(df[0].unique())}
    variant_to_index = {variant: idx for idx, variant in enumerate(df[1].unique())}

    variant_gene_tensor = torch.tensor([
            [variant_to_index[snp], gene_to_index[gene_id]] for snp, gene_id in zip(df[1], df[0])
    ], dtype = torch.long)
    
    torch.save(variant_to_index, os.path.join(OUTPUT_PATH, f'var_index_chr{chromosome}.pt'))
    torch.save(gene_to_index, os.path.join(OUTPUT_PATH, f'gene_index_chr{chromosome}.pt'))
    torch.save(variant_gene_tensor, os.path.join(OUTPUT_PATH, f'variant_gene_map_chr{chromosome}.pt'))
    annot = torch.tensor(np.vstack(annot), dtype = torch.float32)
    torch.save(annot, os.path.join(OUTPUT_PATH, f'annotations_chr{chromosome}.pt'))
    genot = torch.tensor(np.vstack(genotypes), dtype = torch.float32)
    torch.save(genot, os.path.join(OUTPUT_PATH, f'genotypes_chr{chromosome}.pt'))
    return


def main():
    if variant_type == "noncoding":
        cell_df = load_ABC()
    extract = pd.read_csv(VARIANTS_EXTRACT, sep = "\t", header = None)
    extract.columns = ['chromosome','start','end','gene']
    wgsa = pd.read_csv(WGSA_PATH)
    annot, genotypes, variants, genes = load_genes(extract, wgsa)
    save_outputs(annot, genotypes, variants, genes)
    return



if __name__ == "__main__":
    main()
    
    
