################################################################################################
############################ This script loads data for parmigiano ################################
################################################################################################

import utils
import os, sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.decomposition import NMF
import torch
import random

###################################################################################################

def read_data(params):
    '''
    Load phenotype, covariates, and genes
    INPUT: 
        - params: input yaml file loaded as a params dictionary
    OUTPUT: 
        - X: covariate [individual x covariates] dataframe
        - Y: phenotype array
        - genes: dictionary of {chr: [gene1, gene2,...]} to be tested
        - skip: index of individuals to exclude from analysis (if focusing on specific ancestry)
    '''
    skip = None
    XY = pd.read_csv(params['phenotype'])
    if params['ancestry'] != "ALL":
        skip = XY[XY['predicted_ancestry']!=params['ancestry']].index
        XY = XY[XY['predicted_ancestry']==params['ancestry']]
    Y = XY['Diagnosis']
    X = XY[params['covariates']]
    if params['covariate_interactions']:
        X['age_age'] = X['Age'] * X['Age']
        X['age_sex'] = X['Age'] * X['Sex']
        X['age_sex2'] = X['Age'] * X['Sex'] * X['Sex']
    X = utils.scale(X)
    X['intercept'] = 1 # ADDING INTERCEPT TO COVARIATES
    genes = {}
    if "genes" not in params:
        genes = {}
    elif params['genes'] == "all":
        for chr_ in range(1, 23):
            genes[chr_] = list(pd.read_csv(os.path.join(params['gene_matrices'], params['cell'], 'chr' + str(chr_), 'genes.txt'), header = None)[0])
    elif params['genes'] == 'advp':
        advp = pd.read_csv("/gpfs/commons/home/adas/Rare-Variant-AD/data/GWAS_prev/advp_genes.csv")
        genes = advp.groupby('chr')['gene'].agg(list).to_dict()
    elif params['genes'] == 'chr22':
        genes[22] = list(pd.read_csv(os.path.join(params['gene_matrices'], params['cell'], 'chr22', 'genes.txt'), header = None)[0])
        genes[22] = genes[22][0:params['ngenes']]
    elif type(params['genes']) == int:
        genes[params['genes']] = list(pd.read_csv(os.path.join(params['gene_matrices'], params['cell'], 'chr' + str(params['genes']), 'genes.txt'), header = None)[0])
    elif params['genes'] == "skat":
        genes = utils.load_skat_raw(params['skat_path'], params['cell'],'D+B', params['skat_threshold'])
    elif params['genes'] == "staar":
        genes = utils.load_staar_raw(params['staar_path'], params['cell'], params['ngenes'])
    elif params['genes'] == "fst":
        genes = utils.load_fst(params['fst_path'], params['cell'], params['threshold'])
    elif params['genes'] == "random":
        chromosome_genes = {}
        for chr_ in range(1, 23):
            chromosome_genes[chr_] = list(pd.read_csv(os.path.join(params['gene_matrices'], params['cell'], 'chr' + str(chr_), 'genes.txt'), header = None)[0])
        all_genes = [(chromosome, gene) for chromosome, genes in chromosome_genes.items() for gene in genes]
        selected_genes = random.sample(all_genes, params['ngenes'])
        genes = {}
        for chromosome, gene in selected_genes:
            if chromosome not in genes:
                genes[chromosome] = []
            genes[chromosome].append(gene)
    else:
        print("Input gene parameter must be 'advp','all', 'chr22', 'skat', or 'chr20s'. Please enter a valid input.")
        return None
    return X, Y, genes, skip


    
def load_gene(gene, CHRO_NB, params):
    '''
    Load genotype and annotation matrix for a given gene
    '''
    geno = pd.read_feather(os.path.join(params['gene_matrices'],params['cell'],'chr'+str(CHRO_NB), gene + "_geno_imputed.binary"))
    anno = pd.read_csv(os.path.join(params['gene_matrices'], params['cell'],'chr'+ str(CHRO_NB), gene + "_anno_nmf.csv" ))
    anno2 = pd.read_feather(os.path.join(params['gene_matrices'],params['cell'],'chr'+str(CHRO_NB), gene + "_anno.binary")).reset_index()
    anno['splice'] = list(anno2[['SpliceAI_DS_DG', 'SpliceAI_DS_DL', 'SpliceAI_DS_AG', 'SpliceAI_DS_AL']].max(axis = 1)) # Don't use NMF for splicing
    #if params['cell'] == "coding":
    #    anno = anno.drop(['roadmap','enhancer'], axis = 1)
    if geno.shape[1] != anno.shape[0]: # problem with gene input
        return None
    else:
        return geno, anno
    
def load_lof(params, CHRO_NB):
    path = f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/VEP/{params['cell']}/chr{CHRO_NB}_lof.txt"
    lof = pd.read_csv(path, sep = "\t", skiprows = 4).reset_index()
    return lof

def load_genes(params, genes):
    '''
    Loads genotype and functional annotation matrices for genes to be included in model
    INPUT: 
        - params: input yaml file loaded as a params dictionary
        - genes: dictionary of {chr: [gene1, gene2,...]} to be tested
    OUTPUT:
        - Zs: functional annotation dataframe [variants x annotations] with gene mapping included
        - Gs: genotype dataframe [individuals x variants] with matched order of variants to Zs
    '''
    Zs = {}
    Gs = {}
    to_remove = []
    for CHRO_NB in genes:
        if (params['cell'] == "coding") & (params['lof_preds'] == True): lof = load_lof(params, CHRO_NB)
        for gene in tqdm(genes[CHRO_NB]):
            try:
                geno, anno = load_gene(gene, CHRO_NB, params)
                if geno.shape[1] != anno.shape[0]:
                    print("problem with gene dims", gene)
                    to_remove.append([CHRO_NB, gene])
                else:
                    if (params['cell'] == "coding") & (params['lof_preds'] == True):
                        anno['lof'] = np.where(anno['variant_id'].str.split("_").str[0].str.split("chr").str[1].isin(list(lof['level_1'])), 1, 0)
                        print("PCT LOF for ", gene, anno['lof'].sum() / len(anno))
                    Zs[gene] = anno
                    Gs[gene] = geno
            except:
                print("problem reading gene: ", gene)
                to_remove.append([CHRO_NB, gene])
                continue
    for CHRO_NB, gene in to_remove:
        genes[CHRO_NB].remove(gene)
    Gs = pd.concat(Gs, axis=1)
    Gs.columns = Gs.columns.droplevel(0)
    Zs = pd.concat(Zs, axis= 0)
    Zs.index = Zs.index.droplevel(0)
    Zs = Zs.fillna(0) # some genes dont have consistent annotations because they are constant within the gene
    
    if (params['enformer_preds']) & (params['cell'] != 'coding'): 
        enformer = pd.read_csv(params['enformer_path'], sep = "\t")
        delta_columns = [col for col in enformer.columns if 'delta' in col]
        enformer[delta_columns] = enformer[delta_columns].abs()
        Zs = Zs.merge(enformer, left_on = 'variant_id', right_on = 'SNP').drop(['CHR','SNP','BP','A1','A2'], axis = 1)
        
    if type(params['annotations']) == list: 
        Zs = Zs[params['annotations']]
    elif params['annotations'] == 'none':
        Zs = Zs[['intercept']]
    elif params['annotations'] == "NMF_PREP":
        Zs = utils.scale(Zs.set_index(['variant_id', 'TargetGene']))
        Zs['intercept'] = 1 
        print(Zs.columns)
    return Gs, Zs, genes


def load_model(path):
    tau = None
    psi = None
    alpha = None
    d = {'tau': pd.DataFrame(), 'psi': pd.DataFrame(), 'alpha': pd.DataFrame()}
    for iteration in os.listdir(path):
        for variable in ['tau', 'psi', 'alpha']:
            try:
                df = pd.read_csv(os.path.join(path, iteration, variable + '.csv'), index_col = 0).T.reset_index()
                df['iteration'] = int(iteration)
                d[variable] = pd.concat((d[variable], df))
            except:
                continue
    if len(d['tau']) != 0:
        tau = torch.tensor(d['tau'].groupby('index', sort = False).mean()['mean'], dtype = torch.float32)
    else:
        print("Issue loading tau")
    if len(d['psi']) != 0:
        psi = torch.tensor(d['psi'].groupby('index', sort = False).mean()['mean'], dtype = torch.float32)
    else:
        print("Issue loading psi")
    if len(d['alpha']) != 0:
        alpha = torch.tensor(d['alpha'].groupby('index', sort = False).mean()['mean'], dtype = torch.float32)
    else:
        print("Issue loading alpha")
    return tau, psi, alpha
