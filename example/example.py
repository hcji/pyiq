# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 10:23:34 2020

@author: hcji
"""

import pandas as pd
from pyiq.preprocess import preprocess, create_protein_list
from pyiq.iq import create_protein_table


quant_table = pd.read_csv('data/spikeins.csv')
norm_data = preprocess(quant_table)
protein_list = create_protein_list(norm_data)
protein_table = create_protein_table(protein_list)