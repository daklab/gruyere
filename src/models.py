
import pyro
import pyro.distributions as dist
from pyro.nn import PyroSample, PyroModule
import torch


class gruyere(PyroModule):
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
        prior = torch.ones(data.num_anno) / (data.num_anno)
        tau = pyro.sample('tau', dist.Dirichlet(prior))
        if simulate:
            data.wg = w_g
            data.tau = tau
            data.Y['train'] = {}
            data.Y['test'] = {}
        for gene in range(data.num_genes): 
            alpha = pyro.sample(f'alpha_{data.genes[gene]}', dist.Normal(0,1).expand([data.num_cov]).to_event(1))
            beta = ((data.Zs[data.gene_indices==gene].T * data.maf_weights[data.gene_indices == gene]).T.matmul(tau)) * w_g[gene] 
            Gbeta = data.G['train'][:,data.gene_indices==gene].matmul(beta)
            mean = torch.sigmoid(torch.matmul(data.X['train'], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(num_indivs)
            if simulate:
                data.Y['train'][data.genes[gene]] = pyro.sample('obs', dist.Bernoulli(mean))
                Gbeta_test = data.G['test'][:,data.gene_indices==gene].matmul(beta)
                mean_test = torch.sigmoid(torch.matmul(data.X['test'], alpha).reshape(-1,1) + Gbeta_test.reshape(-1,1)).view(len(data.G['test']))
                data.Y['test'][data.genes[gene]] = pyro.sample('obs_test', dist.Bernoulli(mean_test))
            # To Do: Add parameter for non-binary traits -- then pyro.sample from Normal(mean, 1)
            with pyro.plate(f'data_{data.genes[gene]}', num_indivs):
                if params['simulate']:
                    obs = pyro.sample(f'obs_{data.genes[gene]}', dist.Bernoulli(mean), obs = data.Y['train'][data.genes[gene]])
                else:
                    obs = pyro.sample(f'obs_{data.genes[gene]}', dist.Bernoulli(mean), obs = data.Y['train'])
        return data
