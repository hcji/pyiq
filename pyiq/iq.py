# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 10:28:42 2020

@author: hcji
"""

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from tqdm import tqdm


def create_protein_table(protein_list, method = "maxLFQ", N = 3, aggregation_in_log_space = True):
    samples = protein_list[list(protein_list.keys())[0]].columns
    tab = pd.DataFrame(np.full((len(protein_list), len(samples)), np.nan), columns = samples)
    tab.index = protein_list.keys()
    annotations = []
    
    for prot, X in tqdm(protein_list.items()):
        if method == 'maxLFQ':
            try:
                estimate, annotation = maxLFQ(X)
            except:
                estimate, annotation = topN(X, N = N, aggregation_in_log_space = aggregation_in_log_space)
        elif method == 'topN':
            estimate, annotation = topN(X, N = N, aggregation_in_log_space = aggregation_in_log_space)
        elif method == 'meanInt':
            estimate, annotation = meanInt(X, aggregation_in_log_space = aggregation_in_log_space)
        else:
            raise IOError ('invalid method')
        tab.loc[prot, :] = np.array(estimate)
        annotations.append(annotation)
    return tab, annotations


def maxLFQ(X):
    if len(X) == 1:
        estimate = X
        annotation = 'NA'
        return estimate, annotation
    elif np.isnan(X).all(axis=None):
        estimate = None
        annotation = ''
        return estimate, annotation
    else:
        pass
    
    X = np.array(X)
    N = X.shape[1]
    cc = 0
    g = np.full(N, np.nan)
    
    def spread(i):
        g[i] = cc
        for r in range(X.shape[0]):
            if ~np.isnan(X[r, i]):
                for k in range(X.shape[1]):
                    if (~np.isnan(X[r, k])) and (np.isnan(g[k])):
                        spread(k)
    
    def maxLFQ_do(X):
        N = X.shape[1]
        AtA = np.zeros((N, N))
        Atb = np.zeros(N)
        
        for i in range(N-1):
            for j in range(i+1, N):
                r_i_j = np.nanmedian(-X[:,i] + X[:,j])               
                if np.isnan( r_i_j):
                    continue
                
                AtA[i, j] = -1
                AtA[j, i] = -1
                AtA[i, i] = AtA[i, i] + 1
                AtA[j, j] = AtA[j, j] + 1

                Atb[i] = Atb[i] - r_i_j
                Atb[j] = Atb[j] + r_i_j
        
        AA = np.hstack(( np.vstack(((2 * AtA), np.ones(N))), np.expand_dims(np.append(np.ones(N), 0), 1) ))
        bb = np.append(2 * Atb, np.nanmean(X) * N)
        
        estimate, residual = nnls(AA, bb)
        return estimate[range(N)]

    for i in range(N):
        if np.isnan(g[i]):
            cc += 1
            spread(i)
    
    w = np.full(N, np.nan)
    for i in range(1, cc+1):
        ind = np.where(g == i)[0]
        if sum(ind) == 0:
            w[ind] = np.nanmedian(X[:, ind])
        else:
            w[ind] = maxLFQ_do(X[:, ind])
    
    if np.isnan(w).all():
        estimate = w
        annotation = "NA"
    else:
        quantified_samples = np.where(~np.isnan(w))[0]
        if (g[quantified_samples] == g[quantified_samples[0]]).all():
            estimate = w
            annotation = ""
        else:
            estimate = w
            annotation = g
    return estimate, annotation
        

def topN(X, N = 3, aggregation_in_log_space = True):
    if len(X) == 1:
        estimate = X
        annotation = 'NA'
        return estimate, annotation
    
    if aggregation_in_log_space:
        v = np.nanmean(X, axis = 1)
        v_sorted = np.argsort(-v)
        out = np.nanmean(X.iloc[v_sorted[range(min(N, len(v)))],:], axis = 0)
    else:
        XX = 2 ** X
        v = np.nanmean(XX, axis = 1)
        v_sorted = np.argsort(-v)
        out = np.log2(np.nanmean(XX.iloc[v_sorted[range(min(N, len(v)))],:], axis = 0))
    estimate = out
    annotation = ''
    return estimate, annotation


def meanInt(X, aggregation_in_log_space = True):
    if len(X) == 1:
        estimate = X
        annotation = 'NA'
        return estimate, annotation    

    if aggregation_in_log_space:
        out = np.nanmean(X, axis = 0)
    else:
        out = np.log2(np.nanmean(2 ** X, axis = 0))
    estimate = out
    annotation = ''
    return estimate, annotation
