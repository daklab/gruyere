import data_class
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


class SKAT(PyroModule):
    def __init__(self):
        super().__init__()
        
    def forward(self, data, params, simulate = False):
        num_indivs = data.G['train'].shape[0]
        w_g = pyro.sample("w_g", dist.Normal(data.wg_prior,1))
        alpha = pyro.sample('alpha', dist.Normal(0,1).expand([data.num_cov]).to_event(1))   
        lambda_ = ((data.Z.T * data.maf_weights).T.matmul(data.tau)) * w_g
        Z_norm = pyro.sample("Z", dist.Normal(0,1).expand([len(data.maf_weights)]).to_event(1))
        beta = lambda_ * Z_norm
        if simulate:
            data.wg = w_g
            data.alpha = alpha
            data.Z_norm = Z_norm        
        Gbeta = (data.G['train']).matmul(beta) 
        mean = torch.sigmoid(torch.matmul(data.X['train'], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(num_indivs) 
        with pyro.plate('data', num_indivs):
            if simulate:
                data.AD_status['train'] = pyro.sample("obs", dist.Bernoulli(mean))
                return data
            else:
                obs = pyro.sample('obs', dist.Bernoulli(mean), obs = data.AD_status['train'])
        return 
    
class Burden(PyroModule):
    def __init__(self):
        super().__init__()
        
    def forward(self, data, params, simulate = False):
        num_indivs = data.G['train'].shape[0]
        w_g = pyro.sample("w_g", dist.Normal(data.wg_prior,1))
        alpha = pyro.sample('alpha', dist.Normal(0,1).expand([data.num_cov]).to_event(1))   
        lambda_ = ((data.Z.T * data.maf_weights).T.matmul(data.tau)) * w_g
        if simulate:
            data.wg = w_g
            data.alpha = alpha
        Gbeta = (data.G['train']).matmul(lambda_) 
        mean = torch.sigmoid(torch.matmul(data.X['train'], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(num_indivs) 
        with pyro.plate('data', num_indivs):
            if simulate:
                data.AD_status['train'] = pyro.sample("obs", dist.Bernoulli(mean))
                return data
            else:
                obs = pyro.sample('obs', dist.Bernoulli(mean), obs = data.AD_status['train'])
        return 
    
class parmigiano_pergene(PyroModule):
    def __init__(self):
        super().__init__()
        
    def forward(self, data, params, simulate = False):
        num_indivs = data.G['train'].shape[0]
        rho = pyro.sample("rho", dist.Beta(1,0.5))
        w_g = pyro.sample("w_g", dist.Normal(data.wg_prior,1))
        alpha = pyro.sample('alpha', dist.Normal(0,1).expand([data.num_cov]).to_event(1))   
        lambda_ = ((data.Z.T * data.maf_weights).T.matmul(data.tau)) * w_g
        Z_norm = pyro.sample("Z", dist.Normal(0,1).expand([len(data.maf_weights)]).to_event(1))
        beta = rho * lambda_ + (1-rho) * lambda_ * Z_norm
        if simulate:
            data.wg = w_g
            data.alpha = alpha
            data.Z_norm = Z_norm  
            data.rho = rho
        Gbeta = (data.G['train']).matmul(beta) 
        mean = torch.sigmoid(torch.matmul(data.X['train'], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(num_indivs) 
        with pyro.plate('data', num_indivs):
            if simulate:
                data.AD_status['train'] = pyro.sample('obs', dist.Bernoulli(mean))
                print("Data Simulated")
                return data
            else:
                obs = pyro.sample('obs', dist.Bernoulli(mean), obs = data.AD_status['train'])
        return 
    
    
    
class parmigiano(PyroModule):
    def __init__(self):
        super().__init__()
    
    def forward(self, data, params, simulate = False):
        num_indivs = data.G['train'].shape[0]
        w_g = pyro.sample('w_g', dist.Normal(0,1).expand([data.num_genes]).to_event(1))
        rho_g = pyro.sample("rho_g", dist.Beta(0.5,0.5).expand([data.num_genes]).to_event(1))
        #tau_param = pyro.param("tau_param", data.tau)
        #tau = pyro.sample("tau", dist.Uniform(0,1).expand([data.num_anno]).to_event(1))
        prior = torch.ones(data.num_anno) / (data.num_anno)
        tau = pyro.sample('tau', dist.Dirichlet(prior))
        if not params['alpha_gene']:
            alpha = pyro.sample('alpha', dist.Normal(0,1).expand([data.num_cov]).to_event(1))
        if simulate:
            data.wg = w_g
            data.rho = rho_g
            data.tau = tau
            data.AD_status['train'] = {}
        for gene in range(data.num_genes): 
            if params['alpha_gene']:
                alpha = pyro.sample(f'alpha_{data.genes[gene]}', dist.Normal(0,1).expand([data.num_cov]).to_event(1))
            lambda_ = ((data.Zs[data.gene_indices==gene].T * data.maf_weights[data.gene_indices == gene]).T.matmul(tau)) * w_g[gene] 
            Z_norm = pyro.sample(f"Z_{data.genes[gene]}", dist.Normal(0,1).expand([sum(data.gene_indices == gene)]).to_event(1))
            beta = rho_g[gene] * lambda_ + (1-rho_g[gene]) * lambda_ * Z_norm
            Gbeta = data.G['train'][:,data.gene_indices==gene].matmul(beta)
            mean = torch.sigmoid(torch.matmul(data.X['train'], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(num_indivs)
            with pyro.plate(f'data_{data.genes[gene]}', num_indivs):
                if simulate:
                    data.AD_status['train'][data.genes[gene]] = pyro.sample('obs', dist.Bernoulli(mean))
                else:
                    if params['simulate']:
                        obs = pyro.sample(f'obs_{data.genes[gene]}', dist.Bernoulli(mean), obs = data.AD_status['train'][data.genes[gene]])
                    else:
                        obs = pyro.sample(f'obs_{data.genes[gene]}', dist.Bernoulli(mean), obs = data.AD_status['train'])
            with pyro.plate(f'data_{data.genes[gene]}', data.G['test'].shape[0]):
                if simulate:
                    Gbeta = data.G['test'][:,data.gene_indices==gene].matmul(beta)
                    mean = torch.sigmoid(torch.matmul(data.X['test'], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(data.G['test'].shape[0])
                    data.AD_status['test'][data.genes[gene]] = pyro.sample('obs', dist.Bernoulli(mean))
        return data
    
    
    
class parmigiano_separated(PyroModule):
    def __init__(self):
        super().__init__()
    
    def forward(self, data, params, simulate = False):
        num_indivs = data.G['train'].shape[0]
        w_g = pyro.sample('w_g', dist.Normal(0,1).expand([data.num_genes]).to_event(1))
        rho_g = pyro.sample("rho_g", dist.Beta(0.5,0.5).expand([data.num_genes]).to_event(1))
        #tau_param = pyro.param("tau_param", data.tau)
        #tau = pyro.sample("tau", dist.Uniform(0,1).expand([data.num_anno]).to_event(1))
        prior = torch.ones(data.num_anno) / (data.num_anno)
        tau = pyro.sample('tau', dist.Dirichlet(prior))
        psi = pyro.sample('psi', dist.Dirichlet(prior))
        if not params['alpha_gene']:
            alpha = pyro.sample('alpha', dist.Normal(0,1).expand([data.num_cov]).to_event(1))
        if simulate:
            data.wg = w_g
            data.rho = rho_g
            data.tau = tau
            data.psi = psi
            data.AD_status['train'] = {}
        for gene in range(data.num_genes): 
            if params['alpha_gene']:
                alpha = pyro.sample(f'alpha_{data.genes[gene]}', dist.Normal(0,1).expand([data.num_cov]).to_event(1))
            lambda_mu = ((data.Zs[data.gene_indices==gene].T * data.maf_weights[data.gene_indices == gene]).T.matmul(tau)) * w_g[gene] 
            lambda_sigma = ((data.Zs[data.gene_indices==gene].T * data.maf_weights[data.gene_indices == gene]).T.matmul(psi)) * w_g[gene] 
            Z_norm = pyro.sample(f"Z_{data.genes[gene]}", dist.Normal(0,1).expand([sum(data.gene_indices == gene)]).to_event(1))
            beta = rho_g[gene] * lambda_mu + (1-rho_g[gene]) * lambda_sigma * Z_norm
            Gbeta = data.G['train'][:,data.gene_indices==gene].matmul(beta)
            mean = torch.sigmoid(torch.matmul(data.X['train'], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(num_indivs)
            with pyro.plate(f'data_{data.genes[gene]}', num_indivs):
                if simulate:
                    data.AD_status['train'][data.genes[gene]] = pyro.sample('obs', dist.Bernoulli(mean))
                else:
                    if params['simulate']:
                        obs = pyro.sample(f'obs_{data.genes[gene]}', dist.Bernoulli(mean), obs = data.AD_status['train'][data.genes[gene]])
                    else:
                        obs = pyro.sample(f'obs_{data.genes[gene]}', dist.Bernoulli(mean), obs = data.AD_status['train'])
        return data