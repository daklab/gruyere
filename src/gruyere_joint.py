###### This script runs genome-wide gruyere ######

import torch
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
# Loading other scripts
import utils
import load_data
import data_class
import models
import performance



def fit(data, params):
    '''
    Fit gruyere model jointly
    '''
    if params['simulate']:
        data = models.gruyere().forward(data, params, simulate = True)
    model = models.gruyere()
    guide = AutoGuideList(model)
    to_optimize = [] # can have "tau" here if want point estimates (Delta guide)
    guide.add(AutoNormal(poutine.block(model, hide = to_optimize)))
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
        print(simulation_results)
    return posterior_stats, losses

def write_outputs(posterior_stats, losses, params, data, annotations, covariates, train_perf, test_perf):  
    if not os.path.exists(params['output']):
        os.mkdir(params['output'])
    params['output'] = os.path.join(params['output'], 'joint_model')
    if not os.path.exists(params['output']):
        os.mkdir(params['output'])
    np.savetxt(os.path.join(params['output'], "losses.txt"), losses)
    df = pd.DataFrame.from_dict(posterior_stats['tau'], orient = 'index')
    df.columns = annotations
    df.T.to_csv(os.path.join(params['output'], "tau.csv"))
    
    alphas = {key.split('_')[1]: value['mean'] for key, value in posterior_stats.items() if key.startswith('alpha_')}
    alphas = pd.DataFrame(alphas)
    alphas.index = covariates
    alphas.to_csv(os.path.join(params['output'], 'alpha.csv')) 
    
    df = pd.DataFrame.from_dict(posterior_stats['w_g'], orient = 'index')
    df.columns = data.genes
    df = df.T
    df.columns = ['wg', 'wg_std']
    df.to_csv(os.path.join(params['output'], "wg.csv"))

    pd.DataFrame(train_perf).to_csv(os.path.join(params['output'], 'train_performance.csv'))
    try:
        pd.DataFrame(test_perf).to_csv(os.path.join(params['output'], 'test_performance.csv'))
    except: 
        print("No held-out test set")
    return
    
    
def run_gruyere(params):
    X, Y, genes = load_data.read_data(params)
    Gs, Zs = load_data.load_genes(params, genes)
    data = data_class.GenomeWideAD.from_pandas(Gs, Zs, X, Y, params)
    annotations = Zs.columns
    covariates = X.columns
    posterior_stats, losses = fit(data, params)
    train_perf = performance.predict_joint(data, params, posterior_stats, 'train')
    if params['test_prop'] != 0:
        test_perf = performance.predict_joint(data, params, posterior_stats, 'test')
    else: test_perf = None
    write_outputs(posterior_stats, losses, params, data, annotations, covariates, train_perf, test_perf)
    return
    
    
if __name__ == "__main__":
    # Load YAML input arguments
    params_file = sys.argv[1]
    with open(params_file, 'r') as stream:
        params = yaml.safe_load(stream)  
    run_gruyere(params)