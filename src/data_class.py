import pandas as pd
import numpy as np
import os, sys
import torch
import utils
from dataclasses import dataclass, field
from sklearn.model_selection import train_test_split

RANDOM_STATE = 42
    
@dataclass
class GenomeWideAD: # to do: not AD specific
    G: dict
    Zs: torch.Tensor
    X: dict
    gene_indices: torch.Tensor
    genes: pd.Index
    Y: dict
    maf_weights: torch.Tensor
    device: torch.device
    num_genes: int = 0
    num_anno: int = 0
    num_cov: int = 0
    
    def __post_init__(self):
        self.num_genes = max(self.gene_indices) + 1
        self.num_anno = self.Zs.shape[1]
        self.num_cov = self.X['train'].shape[1]
        
    @staticmethod
    def from_pandas(Gs, Zs, X, Y, params, device = "cpu"): 
        gene_indices, genes = pd.factorize(Zs.index.get_level_values("TargetGene"))
        if params['test_prop'] == 0: # If no train test split
            X = {'train': torch.tensor(np.array(X),dtype = torch.float, device = device), 'test': None}
            G = {'train': torch.tensor(np.array(Gs),dtype = torch.float, device = device), 'test': None}
            Y = {'train': torch.tensor(np.array(Y),dtype = torch.float, device = device), 'test': None}
        else:
            X_train, X_test = train_test_split(X, test_size = params['test_prop'], random_state=RANDOM_STATE)
            G_train, G_test = train_test_split(Gs, test_size = params['test_prop'], random_state=RANDOM_STATE)
            Y_train, Y_test = train_test_split(Y, test_size = params['test_prop'], random_state=RANDOM_STATE)
            X_train = torch.tensor(np.array(X_train),dtype = torch.float, device = device)
            X_test = torch.tensor(np.array(X_test),dtype = torch.float, device = device)
            G_train = torch.tensor(np.array(G_train),dtype = torch.float, device = device)
            G_test = torch.tensor(np.array(G_test),dtype = torch.float, device = device)
            Y_train = torch.tensor(np.array(Y_train),dtype = torch.float, device = device)
            Y_test = torch.tensor(np.array(Y_test),dtype = torch.float, device = device)
            G = {'train': G_train, 'test': G_test}
            X = {'train': X_train, 'test': X_test}
            Y = {'train': Y_train, 'test': Y_test}
        
        return GenomeWideAD(
            gene_indices = torch.tensor(gene_indices, dtype = torch.long, device = device), 
            genes = genes,
            Zs = torch.tensor(np.array(Zs), dtype = torch.float, device = device),
            maf_weights = utils.get_weights(Gs), 
            G = G,
            X = X,
            Y = Y,
            device = device 
        )
    
    
@dataclass
class PerGeneAD: 
    G: dict
    Z: torch.Tensor
    X: dict
    Y: dict
    maf_weights: torch.Tensor
    device: torch.device
    num_genes: int = 0
    num_anno: int = 0
    num_cov: int = 0
    
    def __post_init__(self):
        self.num_anno = self.Z.shape[1]
        self.num_cov = self.X['train'].shape[1]
        self.num_genes = 1
        
    @staticmethod
    def from_pandas(Gs, Z, X, Y, params, device = "cpu"): 
        if params['test_prop'] == 0: # If no train test split
            X = {'train': torch.tensor(np.array(X),dtype = torch.float, device = device), 'test': None}
            G = {'train': torch.tensor(np.array(Gs),dtype = torch.float, device = device), 'test': None}
            Y = {'train': torch.tensor(np.array(Y),dtype = torch.float, device = device), 'test': None}
        else:
            X_train, X_test = train_test_split(X, test_size = params['test_prop'], random_state=RANDOM_STATE)
            G_train, G_test = train_test_split(Gs, test_size = params['test_prop'], random_state=RANDOM_STATE)
            Y_train, Y_test = train_test_split(Y, test_size = params['test_prop'], random_state=RANDOM_STATE)
            X_train = torch.tensor(np.array(X_train),dtype = torch.float, device = device)
            X_test = torch.tensor(np.array(X_test),dtype = torch.float, device = device)
            G_train = torch.tensor(np.array(G_train),dtype = torch.float, device = device)
            G_test = torch.tensor(np.array(G_test),dtype = torch.float, device = device)
            Y_train = torch.tensor(np.array(Y_train),dtype = torch.float, device = device)
            Y_test = torch.tensor(np.array(Y_test),dtype = torch.float, device = device)
            G = {'train': G_train, 'test': G_test}
            X = {'train': X_train, 'test': X_test}
            Y = {'train': Y_train, 'test': Y_test}
        
        return PerGeneAD(
            Z = torch.tensor(np.array(Z), dtype = torch.float, device = device),
            maf_weights = utils.get_weights(Gs), 
            G = G,
            X = X,
            Y = Y,
            device = device 
        )
    
   
