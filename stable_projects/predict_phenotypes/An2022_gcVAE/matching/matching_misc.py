#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
import pandas as pd


def get_max_threshold_round(macc_combs):
    """
    Get maximum #Tps for each round data

    Args:
        macc_combs (dict): Dictonary for MACC combinations
    """
    max_threshold = 1
    macc_subs = list(macc_combs.keys())
    for sub in macc_subs:
        if macc_combs[sub]['NumTPs'] > max_threshold:
            max_threshold = macc_combs[sub]['NumTPs']
    return max_threshold


def get_max_threshold(csv_path):
    """
    For MACC_x_bin.csv, get max tps for subjects

    Args:
        csv_path (str): Path for MACC_x_bin.csv
    """
    max_threshold = 1
    df = pd.read_csv(csv_path)
    subjects = np.unique(df.RID)
    for sub in subjects:
        sub_mask = (df.RID == sub)
        sub_data = df[sub_mask]
        if sub_data.shape[0] > max_threshold:
            max_threshold = sub_data.shape[0]
    return max_threshold
