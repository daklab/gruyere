import sklearn
from sklearn import metrics 
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
import numpy as np
import pandas as pd
import torch

def predict_pergene(data, params, posterior_stats, group):
    performance = {}
    # To do: Add per-gene predictions
    return performance
    

def predict_joint(data, params, posterior_stats, group):
    num_indivs = data.G[group].shape[0]
    performance = {}
    performance['AUC'] = {}
    performance['ACC'] = {}
    tau = torch.tensor(posterior_stats['tau']['mean'], dtype = torch.float32)
    w_g = torch.tensor(posterior_stats['w_g']['mean'], dtype = torch.float32)
    for gene in range(data.num_genes):
        alpha = torch.tensor(posterior_stats[f'alpha_{data.genes[gene]}']['mean'], dtype = torch.float32)
        beta = ((data.Zs[data.gene_indices==gene].T * data.maf_weights[data.gene_indices == gene]).T.matmul(tau)) * w_g[gene] 
        Gbeta = data.G[group][:,data.gene_indices==gene].matmul(beta)
        preds = torch.sigmoid(torch.matmul(data.X[group], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(num_indivs)
        if params['simulate']:
            fpr, tpr, thresholds = roc_curve(np.array(data.Y[group][data.genes[gene]]), preds.detach().numpy())
            performance['AUC'][data.genes[gene]] = auc(fpr, tpr)
            performance['ACC'][data.genes[gene]] = ((preds > 0.5).float().detach().numpy() == np.array(data.Y[group][data.genes[gene]])).sum() / len(preds)
        else:
            fpr, tpr, thresholds = roc_curve(np.array(data.Y[group]), preds.detach().numpy())
            performance['AUC'][data.genes[gene]] = auc(fpr, tpr)
            performance['ACC'][data.genes[gene]] = ((preds > 0.5).float().detach().numpy() == np.array(data.Y[group])).sum() / len(preds)
    return performance
