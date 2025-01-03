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
