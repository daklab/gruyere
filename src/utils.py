import pandas as pd
import numpy as np
import os, sys
import torch
<<<<<<< HEAD
import statsmodels
import pyro
import pyro.distributions as dist
from tqdm import tqdm


def scale(df): # Min-max
    return (df-df.min())/ (df.max() - df.min())

def scale_std(df):
    return df.apply(lambda x: x / x.std(), axis=0)
    
def impute(G, type_, device = "cpu"):
    '''
    type_ is the type of imputation, can be either "fixed" or "constant"
    Impute missing genotypes using fixed method from SKAT
    '''
    if type_ == "fixed":
        G = G.to_numpy()
        column_means = np.nanmean(G, axis=0) / 2
        nan_indices = np.isnan(G)
        for col in range(G.shape[1]):
            G[nan_indices[:, col], col] = column_means[col] 
    if type_ == "constant":
        G = G.fillna(0.00005)
    return torch.tensor(np.array(G),dtype = torch.float, device = device)


def get_weights(maf_weights, G, Z):
=======
from tqdm import tqdm
import pyro.distributions as dist



def scale(df):
    # function that min-max scales a dataframe 
    return (df-df.min())/ (df.max() - df.min())

def get_weights(G):
>>>>>>> 3b64698c751033ced518527e22663b311a1c4f33
    '''
    Set variant weights based on MAF and from Beta(1,25) distribution
    This should be run before imputation     
    '''
    d = dist.Beta(1,25)
<<<<<<< HEAD
    if maf_weights == 'observed':
        maf = torch.tensor(G.mean(0)/2, dtype = torch.float)
    else: # if maf_weights == 'gnomad'
        maf = np.exp(-(np.array((Z['gnomAD_genomes_POPMAX_AF'])*(23.718 - 0.6931) + 0.6931)))# reverse the 0-1 and log normalization
        maf = torch.tensor(maf, dtype = torch.float) 
=======
    maf = torch.tensor(G.mean(0)/2, dtype = torch.float) 
>>>>>>> 3b64698c751033ced518527e22663b311a1c4f33
    maf[maf>0.5] = 1 - maf[maf>0.5] # this shouldn't change anything, just checking that AF is the correct direction
    maf[maf>0.05] = 0.05 # in case of any leaks
    weights = torch.exp(d.log_prob(maf)) 
    return weights
<<<<<<< HEAD

def load_skat_raw(path, cell, type_test, threshold):
    skat_df = pd.DataFrame()
    iteration = os.listdir(os.path.join(path, cell))[0]
    for chr_ in os.listdir(os.path.join(path, cell, iteration)):
        try:
            file = os.path.join(path, cell, iteration, chr_)
            df = pd.read_csv(file, sep = "\t")
            df['cell'] = cell
            df['chr'] = chr_.split(".")[0]
            skat_df = pd.concat((skat_df, df))
        except:
            print("Issue reading", file)
    skat_df = skat_df.melt(['Gene', 'cell','chr'])
    skat_df['-log(p)'] = -np.log10(skat_df['value'])
    skat_df['type_test'] = skat_df['variable'].map({'SKAT': 'Dispersion', 'SKATO': 'D+B', 'Burden': 'Burden'})
    skat_df = skat_df[skat_df['type_test']==type_test]
    skat_df = skat_df[skat_df['value']<=threshold]
    skat_df['chr_'] = (skat_df['chr'].str.split("chr").str[1]).astype(int)
    grouped = skat_df.groupby('chr_')['Gene'].apply(list).reset_index()
    genes = {row['chr_']: row['Gene'] for _, row in grouped.iterrows()}
    return genes


def load_staar_raw(path, cell, ngenes):
    '''
    Load the most significant genes from STAAR for analysis 
    Inputs:
        - path: location of STAAR outputs
        - cell: cell to focus on for STAAR results
        - ngenes: number of top genes to analyze
    Outputs:
        - genes: dictionary of chromosome: list of genes that are STAAR significant
    '''
    df = pd.DataFrame()
    for chro in tqdm(os.listdir(os.path.join(path, cell, '1'))):
        chro_df = pd.read_csv(os.path.join(path, cell, '1', chro), sep = "\t")
        chro_df['chr']= chro.split("_")[2].split(".")[0][3:]
        df = pd.concat((df, chro_df))
        df_final = df.sort_values("STAAR-O").head(int(ngenes/2))
        df_final = pd.concat((df_final, df.sort_values("STAAR-O").tail(int(ngenes/2)))) # add control genes, HALF AND HALF right now
        grouped = df_final.groupby('chr')['Gene'].apply(list).reset_index()
        genes = {row['chr']: row['Gene'] for _, row in grouped.iterrows()}
    return genes



def load_fst(path, cell, threshold, type_test = "D+B"):
    '''
    Load the most significant genes from FST for analysis
    Inputs:
        - path: location of FST outputs
        - cell: cell to focus on for FST results
        - ngenes: number of top genes to analyze
    Outputs:
        - genes: dictionary of chromosome: list of genes that are FST significant
    '''
    df = pd.DataFrame()
    for iteration in tqdm(os.listdir(os.path.join(path, cell))):
        if int(iteration) > 10: continue
        for file in os.listdir(os.path.join(path, cell, iteration)):
            if file.startswith("by_region"):
                chr_ = pd.read_csv(os.path.join(path, cell, iteration, file), sep = "\t")
                chr_['chromosome'] = file.split("_")[-1].split(".")[0][3:]
                chr_['iteration'] = int(iteration)
                genes = list(chr_[chr_['type_test']=="Burden"].index) 
                chr_ = chr_[chr_['type_test'] == type_test]
                chr_['Gene'] = genes
                df = pd.concat((df, chr_))  
    df_avg = df[['Gene','minP', 'chromosome']].groupby(['Gene','chromosome']).mean()
    df_avg = df_avg[df_avg['minP'] <= threshold].reset_index()
    grouped = df_avg.groupby('chromosome')['Gene'].apply(list).reset_index()
    genes = {row['chromosome']: row['Gene'] for _, row in grouped.iterrows()}
    return genes

def convert_distributions(variable, distribution, counts, data):
    '''
    This function takes the specified distibutions that variables should be drawn from and converts it to code
    Inputs:
        - variable: This is a string of a learnt variable of parmigiano ie. "tau"
        - distribution: This is a list specifying distribution and optionally parameters ie ["Normal", 0, 1] or ["normal"]
        - data: This is a data_class object with genotype, annotation, and phenotype data
    Outputs:
        - A pyro distribution with parameters
    '''
    if distribution[variable][0].lower() == "dirichlet":
        prior = torch.ones(counts[variable]) / (counts[variable])
        return dist.Dirichlet(prior)
    elif distribution[variable][0].lower() == "normal":
        if len(distribution[variable]) == 1:
            return dist.Normal(0.,1.)
        else:
            return dist.Normal(distribution[variable][1], distribution[variable][2])
    elif distribution[variable][0].lower() == "gamma":
        if len(distribution[variable]) == 1:
            return dist.Gamma(1.,2.)
        else:
            return dist.Gamma(distribution[variable][1], distribution[variable][2])
    elif distribution[variable][0].lower() == "uniform":
        if len(distribution[variable]) == 1:
            return dist.Uniform(0.,1.)
        else:
            return dist.Uniform(distribution[variable][1], distribution[variable][2])
    return "Distribution not found"

=======
>>>>>>> 3b64698c751033ced518527e22663b311a1c4f33
