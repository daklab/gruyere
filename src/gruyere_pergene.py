###### This script runs per-gene gruyere ######

import torch
from tqdm import tqdm
import numpy as np
import pandas as pd
import time
import os, sys
from dataclasses import dataclass, field
import yaml
from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm
from scipy import stats
# Loading other scripts
import utils
import load_data
import data_class
import models
import performance

def write_outputs(params, chromosome, wgs, preds):
    if not os.path.exists(params['output']):
        os.mkdir(params['output'])
    params['output'] = os.path.join(params['output'], 'pergene_regression')
    if not os.path.exists(params['output']):
        os.mkdir(params['output'])
    wgs.to_csv(os.path.join(params['output'], "pvals_chr" + str(chromosome) + ".csv"))
    preds.to_csv(os.path.join(params['output'], 'preds_chr' + str(CHRO_NB) + ".csv"))
    return


def run_gruyere(params, chromosome):
    X, Y, _ = load_data.read_data(params)
    Gs, Zs = load_data.load_genes_chromosome(params, chromosome)
    tau = load_data.load_tau(params) # Load globally-fit tau
    # Fit covariate-only model:
    model_reduced = sm.GLM(np.array(Y), np.array(X), family = sm.families.Binomial(link = sm.families.links.logit())).fit()
    wgs = {}
    preds = {}
    for gene in tqdm(Zs.index.get_level_values("TargetGene").unique()):
        G = Gs[Gs.filter(regex=f'^{gene}_').columns]
        Z = Zs[Zs.index.get_level_values("TargetGene")==gene]
        data = data_class.PerGeneAD.from_pandas(G, Z, X, Y, params)
        beta = (data.Z.T * data.maf_weights).T.matmul(tau)
        Gbeta = data.G['train'] @ beta
        X_input = torch.cat((data.X['train'], Gbeta.reshape(-1,1)), 1)
        try:
            model = sm.Logit(np.array(data.Y['train']), np.array(X_input)).fit(maxiter = 200)
            df = model.df_model - model_reduced.df_model
            lr_stat = 2 * (model.llf - model_reduced.llf)
            pval = stats.chi2.sf(lr_stat, df)
            true = model.params[-1]
            preds[gene] = model.predict(np.array(X_input))
            wgs[gene] = {"pval":pval, 'coef':true}
        except:
            print("Issue with logistic regression for gene:", gene)

    wgs = pd.DataFrame.from_dict(wgs).T
    preds = pd.DataFrame(preds)
    # to do: Add auc accuracy for predictions
    write_outputs(params, chromosome, wgs, preds)
    return
    
    
if __name__ == "__main__":
    # Load YAML input arguments
    params_file = sys.argv[1]
    chromosome = sys.argv[2]
    with open(params_file, 'r') as stream:
        params = yaml.safe_load(stream)  
    run_gruyere(params, chromosome)