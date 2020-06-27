# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:37:26 2020

@author: hcji
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm


def preprocess(quant_table, primary_id = "PG.ProteinGroups",
                       secondary_id = ["EG.ModifiedSequence", "FG.Charge", "F.FrgIon", "F.Charge"],
                       sample_id = "R.Condition",
                       intensity_col = "F.PeakArea",
                       median_normalization = True,
                       log2_intensity_cutoff = 0):
    
    second_id = quant_table[secondary_id[0]].astype(str)
    if len(secondary_id) > 1:
        for i in range(1,len(secondary_id)):
            second_id = second_id.str.cat(quant_table[secondary_id[i]].astype(str), sep='_')
    
    d = pd.DataFrame({'protein_list': quant_table[primary_id],
                      'sample_list': quant_table[sample_id],
                      'quant': np.log2(quant_table[intensity_col]),
                      'id': second_id})
    
    # remove nan values
    keep = np.where(~np.isnan(d['quant']))[0]
    d = d.iloc[keep,:]
    
    # remove low intensity features
    a = plt.hist(d['quant'], 100)
    plt.xlabel('log2 intensity')
    plt.vlines(log2_intensity_cutoff, 0, np.max(a[0]), color='red')
    plt.show()
    
    keep = np.where(d['quant'] > log2_intensity_cutoff)[0]
    d = d.iloc[keep,:]
    
    # normalization
    m = []
    dl = []
    for sample, sub in d.groupby('sample_list'):
        dl.append(sub)
        m.append(np.nanmedian(sub['quant']))
    if median_normalization:
        f = np.mean(m) - m
        dl_n = [dl[x]['quant'] + f[x] for x in range(len(m))]
        for x in range(len(dl)):
            dl[x]['quant'] = dl_n[x]
        d = pd.concat(dl)
    
    return d



def create_protein_list(preprocessed_data, missing_value_filter = 0.5):
    samples = np.unique(preprocessed_data['sample_list'])
    p_list = {}
    
    for prot, sub in tqdm(preprocessed_data.groupby('protein_list')):
        idx = np.unique(sub['id'])
        m = pd.DataFrame(np.full((len(idx), len(samples)), np.nan))
        m.index = idx
        m.columns = samples
        for r in sub.index:
            col = sub['sample_list'][r]
            row = sub['id'][r]
            val = sub['quant'][r]
            if np.isnan(m.loc[row, col]):
                m.loc[row, col] = val
            else:
                raise ValueError('Error: duplicate entry')
        missR = [len(np.where(np.isnan(m.loc[s,:]))[0]) for s in m.index]
        missR = np.array(missR) / m.shape[1]
        keep = np.where(missR < missing_value_filter)[0]
        
        if len(keep) < 1:
            continue
        else:
            m = m.iloc[keep,:]    
            p_list[prot] = m
    return p_list