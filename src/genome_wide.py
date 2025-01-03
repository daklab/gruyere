import torch
import torch.nn as nn
import torch.nn.functional as F
import pyro
from pyro import poutine
import pyro.distributions as dist
from pyro.nn import PyroSample, PyroModule
from pyro.infer.autoguide import AutoDiagonalNormal, AutoGuideList, AutoDelta, AutoMultivariateNormal, AutoNormal
from pyro.infer import SVI, Trace_ELBO, RenyiELBO
from pyro.infer import Predictive
from torch.distributions import constraints
from tqdm import tqdm
import numpy as np
import pandas as pd
import time
import os, sys
from dataclasses import dataclass, field
import yaml
import sklearn
from sklearn import metrics 
import json
import pickle
import torch.nn.functional as F
from sklearn.decomposition import NMF
import utils
import data_class
from dataclasses import fields
import load_data
import models
from sklearn.linear_model import LogisticRegression
import performance_metrics


##### INPUT ARGUMENTS ####
params_file = sys.argv[1]
with open(params_file, 'r') as stream:
    params = yaml.safe_load(stream)  
try:
    iteration = sys.argv[2]
except IndexError:
    iteration = False
try:
    annotation_perm = str(sys.argv[3])
except IndexError:
    annotation_perm = False

def fit(data, params):
    if params['simulate']:
        data = getattr(models, params['model'])().forward(data, params, True)
    model = getattr(models, params['model'])()
    guide = AutoGuideList(model)
    to_optimize = ['rho', 'tau']
    guide.add(AutoDiagonalNormal(poutine.block(model, hide = to_optimize)))
    guide.add(AutoDelta(poutine.block(model, expose = to_optimize)))
    adam = pyro.optim.Adam({"lr": params['lr']})
    svi = SVI(model, guide = guide, optim = adam, loss=Trace_ELBO()) 
    pyro.clear_param_store()
    losses = []
    for j in tqdm(range(params['epochs'])): 
        loss = svi.step(data, params)
        losses.append(float(loss))
    guide.requires_grad_(False)
    predictive = Predictive(model, guide=guide, num_samples=10) 
    samples = predictive(data, params)
    posterior_stats = {k:{'mean':np.array(torch.mean(v, 0)[0]),
                         'std': np.array(torch.std(v,0)[0])} for k,v in samples.items() if "obs" not in k}
    print("Finished training")
    if params['simulate']:
        simulation_results = {}
        simulation_results['tau_R'] = np.corrcoef(data.tau, posterior_stats['tau']['mean'])[0,1]
        simulation_results['wg_R'] = np.corrcoef(data.wg, posterior_stats['w_g']['mean'])[0,1]
        simulation_results['rho_R'] = np.corrcoef(data.rho, posterior_stats['rho_g']['mean'])[0,1]
        print("SIMULATION RESULTS", simulation_results)
    return posterior_stats, losses
    
    
def write_outputs(posterior_stats, losses, params, data, annotations, covariates, train_perf, test_perf):
    np.savetxt(os.path.join(params['output_path'], "losses.txt"), losses)
    df = pd.DataFrame.from_dict(posterior_stats['tau'], orient = 'index')
    df.columns = annotations
    df = df.T
    if params['simulate']:
        df['true_tau'] = np.array(data.tau)
    df.to_csv(os.path.join(params['output_path'], "tau.csv"))
    
    try:
        df = pd.DataFrame.from_dict(posterior_stats['psi'], orient = 'index')
        df.columns = annotations
        df = df.T
        if params['simulate']:
            df['true_psi'] = np.array(data.psi)
        df.to_csv(os.path.join(params['output_path'], "psi.csv"))
    except:
        print("No psi")
    
    if not params['alpha_gene']:
        try:
            df = pd.DataFrame.from_dict(posterior_stats['alpha'], orient = "index")
            df.columns = covariates
            df = df.T
            df.to_csv(os.path.join(params['output_path'], "alpha.csv"))
        except:
            print("issue with covariates")
    
    df = pd.DataFrame.from_dict(posterior_stats['w_g'], orient = 'index')
    df.columns = data.genes
    df = df.T
    df.columns = ['wg', 'wg_std']
    if params['simulate']:
        df['true_wg'] = np.array(data.wg)
    df2 = pd.DataFrame.from_dict(posterior_stats['rho_g'], orient = 'index')
    df2.columns = data.genes
    df2 = df2.T
    df2.columns = ['rho_g', 'rho_std']
    if params['simulate']:
        df2['true_rho'] = np.array(data.rho)
    df = pd.concat((df, df2), axis = 1)
    df.to_csv(os.path.join(params['output_path'], "wg_rho.csv"))
    
    #pd.DataFrame(train_perf).to_csv(os.path.join(params['output_path'], 'train_performance.csv'))
    #pd.DataFrame(test_perf).to_csv(os.path.join(params['output_path'], 'test_performance.csv'))
    return
    
    
    
def main():
    print(params)
    if params['simulate']: 
        params['output_path'] = os.path.join(params['output_path'], 'simulation', params['output'])
    elif params['permute_anno']: 
        params['output_path'] = os.path.join(params['output_path'], 'permutation', params['output'])
        if not os.path.exists(params['output_path']):
            os.mkdir(params['output_path'])
        params['output_path'] = os.path.join(params['output_path'], annotation_perm)
    else:
        params['output_path'] = os.path.join(params['output_path'], 'true', params['output'])
    if not os.path.exists(params['output_path']):
        os.mkdir(params['output_path'])
    if iteration:
        params['output_path'] = os.path.join(params['output_path'], str(iteration))
        if not os.path.exists(params['output_path']):
            os.mkdir(params['output_path'])
    X, Y, genes, skip = load_data.read_data(params)
    Gs, Zs, genes = load_data.load_genes(params, genes)
    if (params['permute_anno']) & (type(annotation_perm)==str):
        Zs[annotation_perm] = 1 - np.random.permutation(Zs[annotation_perm].values)
    data = data_class.GenomeWideAD.from_pandas(Gs, Zs, X, Y, params)
    annotations = Zs.columns
    covariates = X.columns
    gene_var = Zs.reset_index()[['variant_id','TargetGene']]
    posterior_stats, losses = fit(data, params)
    train_perf = performance_metrics.predict_joint(data, params, posterior_stats, 'train')
    test_perf = performance_metrics.predict_joint(data, params, posterior_stats, 'test')
    train_perf = None
    test_perf = None
    write_outputs(posterior_stats, losses, params, data, annotations, covariates, train_perf, test_perf)
    return

if __name__ == "__main__":
    main()