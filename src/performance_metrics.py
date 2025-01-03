import sklearn
from sklearn import metrics 
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc
import numpy as np
import pandas as pd
import torch


def calculate_FDR(threshold, preds, data, group):
    predicted_positive = (preds.detach().numpy() > threshold)
    actual_positive = np.array(data.AD_status[group]) == 1
    # True Positives (TP) and False Positives (FP)
    TP = ((predicted_positive) & (actual_positive)).sum()
    FP = ((predicted_positive) & (~actual_positive)).sum()
    if TP + FP > 0:
        FDR = FP / (TP + FP)
    else:
        FDR = 0  # Handle division by zero case
    return FDR

def compute_confusion_metrics(threshold, preds, data, group):
    y_true = np.array(data.AD_status[group])
    y_pred = preds.detach().numpy() > threshold
    TP = ((y_pred) & (y_true == 1)).sum()
    FP = ((y_pred) & (y_true == 0)).sum()
    FN = ((~y_pred) & (y_true == 1)).sum()
    TN = ((~y_pred) & (y_true == 0)).sum()
    return TP, FP, FN, TN


def predict_pergene(data, params, posterior_stats, group):
    performance = {}
    num_indivs = data.G[group].shape[0]
    w_g = torch.tensor(posterior_stats['w_g']['mean'], dtype = torch.float32)
    rho_g = torch.tensor(posterior_stats['rho']['mean'], dtype = torch.float32)
    Z_norm = torch.tensor(posterior_stats["Z"]['mean'], dtype = torch.float32)
    alpha = torch.tensor(posterior_stats['alpha']['mean'], dtype = torch.float32)
    lambda_ = ((data.Z.T * data.maf_weights).T.matmul(data.tau)) * w_g
    beta = rho_g * lambda_ + (1-rho_g) * lambda_ * Z_norm
    Gbeta = (data.G[group]).matmul(beta) 
    preds = torch.sigmoid(torch.matmul(data.X[group], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(num_indivs) 
    fpr, tpr, thresholds = roc_curve(np.array(data.AD_status[group]), preds.detach().numpy())
    performance['AUC'] = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(np.array(data.AD_status[group]), preds.detach().numpy())
    performance['PRC'] = auc(recall, precision)
    performance['ACC'] = ((preds > 0.5).float().detach().numpy() == np.array(data.AD_status[group])).sum() / len(preds)
    for threshold in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        performance[f'FDR_{threshold}'] = calculate_FDR(threshold, preds, data, group)
    TP, FP, FN, TN = compute_confusion_metrics(0.5, preds, data, group)
    performance['Precision'] = TP / (TP + FP) if TP + FP > 0 else 0
    performance['Recall'] = TP / (TP + FN) if TP + FN > 0 else 0
    performance['Specificity'] = TN / (TN + FP) if TN + FP > 0 else 0
    performance['Balanced Accuracy'] = (performance['Recall'] + performance['Specificity']) / 2
    performance['F1_Score'] = (2 * performance['Precision'] * performance['Recall'] /
                                (performance['Precision'] + performance['Recall'])) if (performance['Precision'] + performance['Recall']) > 0 else 0
    performance['FNR'] = FN / (TP + FN) if TP + FN > 0 else 0
    # Matthews Correlation Coefficient (MCC)
    numerator = (TP * TN) - (FP * FN)
    denominator = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))**0.5
    performance['MCC'] = numerator / denominator if denominator > 0 else 0

    return performance
    

def predict_joint(data, params, posterior_stats, group):
    num_indivs = data.G[group].shape[0]
    performance = {}
    performance['AUC'] = {}
    performance['ACC'] = {}
    tau = torch.tensor(posterior_stats['tau']['mean'], dtype = torch.float32)
    w_g = torch.tensor(posterior_stats['w_g']['mean'], dtype = torch.float32)
    rho_g = torch.tensor(posterior_stats['rho_g']['mean'], dtype = torch.float32)
    for gene in range(data.num_genes):
        Z_norm = torch.tensor(posterior_stats[f"Z_{data.genes[gene]}"]['mean'], dtype = torch.float32)
        if params['alpha_gene']:
            alpha = torch.tensor(posterior_stats[f'alpha_{data.genes[gene]}']['mean'], dtype = torch.float32)
        else:
            alpha = torch.tensor(posterior_stats['alpha']['mean'], dtype = torch.float32)
        lambda_ = ((data.Zs[data.gene_indices==gene].T * data.maf_weights[data.gene_indices == gene]).T.matmul(tau)) * w_g[gene] 
        beta = rho_g[gene]  * lambda_ + (1-rho_g[gene]) * lambda_ * Z_norm
        Gbeta = data.G[group][:,data.gene_indices==gene].matmul(beta)
        preds = torch.sigmoid(torch.matmul(data.X[group], alpha).reshape(-1,1) + Gbeta.reshape(-1,1)).view(num_indivs)
        if params['simulate']:
            print("GROUP", group)
            print("data.AD_status[group] shape:", data.AD_status[group].shape)
            print("data.AD_status[group] type:", type(data.AD_status[group]))
            print("data.genes[gene]:", data.genes[gene])
            print("data.genes[gene] type:", type(data.genes[gene]))

            print("data.AD_status[group][data.genes[gene]]", data.AD_status[group][data.genes[gene]])
            print("SHAPE", np.array(data.AD_status[group][data.genes[gene]]).shape)
            print("SHAPE", preds.detach().numpy().shape)
            print("UNIQUE AD STATUS",np.unique(data.AD_status[group][data.genes[gene]]))
            print("RANGE PREDS",np.min(preds.detach().numpy()), np.max(preds.detach().numpy()))
            print(np.isnan(data.AD_status[group][data.genes[gene]]).any())
            print(np.isnan(preds.detach().numpy()).any())
            
            y_true = np.array(data.AD_status[group][data.genes[gene]])
            y_score = preds.detach().numpy()
            fpr, tpr, thresholds = roc_curve(y_true, y_score)
            performance['AUC'][data.genes[gene]] = auc(fpr, tpr)
            performance['ACC'][data.genes[gene]] = ((preds > 0.5).float().detach().numpy() == np.array(data.AD_status[group][data.genes[gene]])).sum() / len(preds)
        else:
            fpr, tpr, thresholds = roc_curve(np.array(data.AD_status[group]), preds.detach().numpy())
            performance['AUC'][data.genes[gene]] = auc(fpr, tpr)
            performance['ACC'][data.genes[gene]] = ((preds > 0.5).float().detach().numpy() == np.array(data.AD_status[group])).sum() / len(preds)
    return performance
