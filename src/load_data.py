import utils
import os, sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import torch


def read_data(params):
    '''
    Load phenotypes and covariates
    INPUT: 
        - params: input yaml file loaded as a params dictionary
    OUTPUT: 
        - X: covariate [individual x covariates] dataframe
        - Y: phenotype array
    '''
    XY = pd.read_csv(params['XY'])
    X = XY.drop(['Diagnosis'], axis = 1)
    X = utils.scale(X)
    X['intercept'] = 1 # Add intercept to covariates
    Y = XY['Diagnosis']
    with open(params['genes'], "r") as file:
        genes = file.read().splitlines()  
    return X, Y, genes



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
    Zs = pd.DataFrame()
    Gs = pd.DataFrame()
    for chro in sorted(os.listdir(params['G'])):
        if chro.endswith(".csv") or chro.endswith(".txt"):
            chro_df = pd.read_csv(os.path.join(params['G'], chro), index_col = 0)
            Gs = pd.concat((Gs, chro_df.loc[chro_df.index.str.split("_").str[0].isin(genes)]))
            chro_df = pd.read_csv(os.path.join(params['Z'], chro), index_col = 0)
            Zs = pd.concat((Zs, chro_df.loc[chro_df.index.str.split("_").str[0].isin(genes)]))
    Gs = Gs.T
    Zs = Zs.fillna(0) # Can adjust imputing missing annotations accordingly (we have no missingness)
    if ("Intercept" not in Zs.columns) or ("intercept" not in Zs.columns):
        Zs['intercept'] = 1
    Zs['TargetGene'] = Zs.index.str.split("_").str[0]
    Zs.set_index('TargetGene', append=True, inplace=True)
    if Gs.shape[1] != Zs.shape[0]: 
        print("Error loading genes. Dimensions of genotype and annotation variants do not match.")
        return None
    return Gs, Zs

def load_genes_chromosome(params, chromosome):
    '''
    Loads genotype and functional annotation matrices for genes in a chromosome
    INPUT: 
        - params: input yaml file loaded as a params dictionary
        - genes: dictionary of {chr: [gene1, gene2,...]} to be tested
    OUTPUT:
        - Zs: functional annotation dataframe [variants x annotations] with gene mapping included
        - Gs: genotype dataframe [individuals x variants] with matched order of variants to Zs
    '''
    Gs = pd.read_csv(os.path.join(params['G'], 'chr'+ str(chromosome) + '.csv'), index_col = 0)
    Zs = pd.read_csv(os.path.join(params['Z'], 'chr'+ str(chromosome) + '.csv'), index_col = 0)
    Gs = Gs.T
    Zs = Zs.fillna(0)
    if ("Intercept" not in Zs.columns) or ("intercept" not in Zs.columns):
        Zs['intercept'] = 1
    Zs['TargetGene'] = Zs.index.str.split("_").str[0]
    Zs.set_index('TargetGene', append=True, inplace=True)
    if Gs.shape[1] != Zs.shape[0]: 
        print("Error loading genes. Dimensions of genotype and annotation variants do not match.")
        return None
    return Gs, Zs



def load_tau(params):
    try:
        tau = torch.tensor(pd.read_csv(os.path.join(os.path.join(params['output'],'joint_model','tau.csv')), index_col = 0)['mean'], dtype = torch.float32)
        return tau
    except:
        print("Issue loading tau. Please make sure you have trained joint model first.")
        return None
