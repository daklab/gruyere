import pandas as pd
import numpy as np
import os, sys
import torch
from tqdm import tqdm
import pyro.distributions as dist



def scale(df):
    # function that min-max scales a dataframe 
    return (df-df.min())/ (df.max() - df.min())

def get_weights(G):
    '''
    Set variant weights based on MAF and from Beta(1,25) distribution
    This should be run before imputation     
    '''
    d = dist.Beta(1,25)
    maf = torch.tensor(G.mean(0)/2, dtype = torch.float) 
    maf[maf>0.5] = 1 - maf[maf>0.5] # this shouldn't change anything, just checking that AF is the correct direction
    maf[maf>0.05] = 0.05 # in case of any leaks
    weights = torch.exp(d.log_prob(maf)) 
    return weights
