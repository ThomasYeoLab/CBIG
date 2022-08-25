#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
import pandas as pd
from itertools import combinations
from utils.misc import save_pkl


def get_within_dataset_combs(data_path, save_path):
    """
    Generate all combinations within dataset for measures we want to match

    Args:
        data_path (str): Path for input data
        save_path (str): Path for saving output
    """
    # load dataset
    columns = ['RID', 'EXAMDATE', 'AGE', 'SEX', 'MMSE', 'DX']
    df = pd.read_csv(data_path, usecols=columns)
    subjects = np.unique(df.RID)
    comb_dict = dict()
    for sub in sorted(subjects):
        sub_data = df[(df.RID == sub)]
        nb_tps = sub_data.shape[0]
        comb_dict[sub] = {}
        comb_dict[sub]['NumTPs'] = nb_tps
        comb_dict[sub]['SEX'] = sub_data.loc[:, 'SEX'].values[0]
        comb_dict[sub]['TPs'] = []
        dates = np.sort(sub_data['EXAMDATE'].values)
        for date in dates:
            date_mask = (sub_data.EXAMDATE == date)
            tp_list = []
            tp_list.append(date)
            tp_list.append(sub_data.loc[date_mask, 'AGE'].values[0])  # AGE
            tp_list.append(sub_data.loc[date_mask, 'MMSE'].values[0])  # MMSE
            tp_list.append(sub_data.loc[date_mask, 'DX'].values[0])  # DX
            comb_dict[sub]['TPs'].append(tp_list)
        # generate all combinations within dataset
        comb_dict[sub]['AGE'] = gen_comb_matrix(
            comb_dict[sub]['TPs'], nb_tps, index=1)
        comb_dict[sub]['MMSE'] = gen_comb_matrix(
            comb_dict[sub]['TPs'], nb_tps, index=2)
        comb_dict[sub]['DX'] = gen_comb_matrix(
            comb_dict[sub]['TPs'], nb_tps, index=3)
    # save to a .pkl file
    save_pkl(comb_dict, save_path)


def gen_comb_matrix(tps_list, nb_tps, index):
    """
    Generate a combination maxtrix for given parameters

    Args:
        tps_list (list): List for time points
        nb_tps (int): Number of time points
        index (int): Index for Age/MMSE/DX
    """
    sub_measure = dict()
    threshold = min(nb_tps, 4)  # 4 is the maximum visits of AIBL
    for tps in range(threshold, 0, -1):
        row_comb = list(combinations([x for x in range(nb_tps)], tps))
        meas_matrix = np.zeros((len(row_comb), tps))
        # fill the matrix
        for row_id in range(len(row_comb)):
            for tp in range(tps):
                meas_matrix[row_id, tp] =\
                    tps_list[row_comb[row_id][tp]][index]
        sub_measure[tps] = {}
        sub_measure[tps]['Comb'] = row_comb
        sub_measure[tps]['Matrix'] = meas_matrix

    return sub_measure
