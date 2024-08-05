#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np


def data_pred_metric(RID_column, pred_prob, truth, threshold):
    """
    Compute subject level dataset prediction accuracy
    """
    # get prediction
    pred = (pred_prob >= threshold) * 1.
    subs = np.unique(RID_column)
    subs, ind = np.unique(RID_column, return_index=True)
    subs = subs[np.argsort(ind)]
    nb_subs = subs.shape[0]
    # Copy RID_column along axis 1, generate a matrix #Tps by #Subs
    RID_matrix = np.tile(RID_column, (nb_subs, 1)).transpose()
    mask = (~(RID_matrix == subs)) * 1.  # 0 means select, 1 means not select
    # copy pred along axis 1, generate a matrix #Tps by #Subs
    pred_matrix = np.tile(pred, (nb_subs, 1)).transpose()
    # 0 means wrong prediction, 1 means correct prediction
    pred_acc_matrix = ((pred_matrix.transpose() == truth) * 1.0).transpose()
    acc_vector = np.mean(np.ma.masked_array(pred_acc_matrix, mask), axis=0)
    acc_vector = np.array(acc_vector)
    return acc_vector, np.mean(acc_vector)
