#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Written by Valeria Kebets and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""
# This function averages a 419x419 correlation matrix
# between/within networks (18x18)
# Based on function written by Csaba Orban

# Inputs:
# The input is a .mat file containing 419x419 matrix (corrmat_fname),
# which will be averaged within/between 18 networks
# - corrmat_fname : path and name of .mat file containing correlation
# matrix to average
# - corrmat_var : name of variable containing correlation matrix
# (in the input .mat file)
# - working_dir : directory where .csv and .txt files for re-indexing
# parcels are located
# - out_fname : path and name of .mat file created as output

# Outputs:
# The output is a .mat file (out_fname) containing a 18x18 matrix
# - thisLC_RSFC_loadings_avg : name of variable containing
# averaged (18x18x) correlation matrix

import numpy as np
import pandas as pd
import sys
import scipy.io as sio

corrmat_fname = sys.argv[1]
corrmat_var = sys.argv[2]
working_dir = sys.argv[3]
out_fname = sys.argv[4]

# load index & labels
newindex = pd.read_csv(
    working_dir + '/' + 'CBIG_Schaefer419_index.csv', header=None).values.T[0]

labels = list(
    pd.read_csv(
        working_dir + '/' + 'networks_label_sorted_NEWORDERING.txt',
        header=None,
        index_col=0).index)
labels[-19:] = ['SC_' + x for x in labels[-19:]]

# networks names
cortex_keywords = [
    'TempPar', 'DefaultC', 'DefaultB', 'DefaultA', 'ContC', 'ContB', 'ContA',
    'Limbic_TempPole', 'Limbic_OFC', 'SalVentAttnB', 'SalVentAttnA',
    'DorsAttnB', 'DorsAttnA', 'SomMotB', 'SomMotA', 'VisPeri', 'VisCent', 'SC'
]

# load 3D matrix
CM3D = np.nan_to_num(sio.loadmat(corrmat_fname)[corrmat_var])
dim3 = CM3D.shape[2]  # number of components/samples

# create 3D matrix where to move averaged CMs from all components/samples
CM = np.zeros((18, 18, dim3))

# for loop around components
for i in np.arange(0, dim3):
    corrmat = CM3D[:, :, i]
    corrmat = pd.DataFrame(corrmat)

    # reindex
    corrmat_len = corrmat.shape[0]
    corrmat = corrmat.reindex_axis(newindex, axis=0)
    corrmat = corrmat.reindex_axis(newindex, axis=1)

    corrmat.index = labels
    corrmat.columns = labels

    corrmat2 = corrmat.copy()  # contains averaged FC in 419x419 matrix
    corrmat3 = pd.DataFrame(
        0, index=np.arange(0, 18),
        columns=np.arange(0, 18))  # contains averaged FC in 18x18 matrix

    for network_1st in cortex_keywords:
        for network_2nd in cortex_keywords:
            network_labels_1st = [x for x in labels if network_1st in x]
            network_labels_2nd = [x for x in labels if network_2nd in x]
            if network_labels_1st == network_labels_2nd:
                # diagonals have to be masked before averaging
                RSFC_mask = np.triu(
                    np.ones((len(network_labels_1st),
                             len(network_labels_2nd))),
                    k=1).astype(np.bool)

                corrmat3.loc[cortex_keywords.index(network_1st),
                             cortex_keywords.index(network_2nd)] = np.nanmean(
                                 corrmat2.
                                 loc[network_labels_1st, network_labels_2nd].
                                 values[RSFC_mask])

            else:
                corrmat3.loc[cortex_keywords.index(network_1st),
                             cortex_keywords.index(network_2nd)] = np.nanmean(
                                 corrmat2.
                                 loc[network_labels_1st, network_labels_2nd])

    CM[:, :, i] = corrmat3  # move matrix

# save to .mat
sio.savemat(out_fname, {'thisLC_RSFC_loadings_avg': CM})
