#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import numpy as np
import scipy.io as sio
from config import config


def get_phes(folder):
    '''get the phenotypes of test meta-set

    Args:
        folder (str): folder path that stores data input files

    Returns:
        List: list of test meta-set phenotypes
    '''

    txt = os.path.join(folder, 'HCP_diff_roi_final_phe_list.txt')

    phes = np.genfromtxt(txt, dtype=str)
    print(phes)
    return phes


def get_haufe_pred():
    '''get data for PNF calculation for KRR ground truth and 100 shot

    Args:
        None

    Returns:
        None
    '''

    base_dir = config.REP_DIR

    haufe_dir = os.path.join(base_dir, 'output_KRR_classical_haufe_100')
    pred_haufe_dict = {}
    all_train_dir = os.path.join(base_dir, 'output_KRR_classical_haufe_all')
    pred_all_dict = {}

    output_dir = os.path.join(base_dir, 'output_intermediate')
    os.makedirs(output_dir, exist_ok=True)

    input_dir = os.path.join(os.environ['CBIG_REPDATA_DIR'],
                             'stable_projects/predict_phenotypes/He2022_MM')
    phes_tes = get_phes(input_dir)
    rngs = 100
    ks = [10, 20, 50, 100, 200]

    for _, phe in enumerate(phes_tes):
        pred_haufe_dict[phe] = np.zeros((rngs, 100))
        pred_all_dict[phe] = np.zeros((100))

    # get split file for k shot
    for _, phe in enumerate(phes_tes):
        for ik, k in enumerate(ks):
            for rng in range(1, rngs + 1):
                if k == 100:
                    folder = os.path.join(
                        haufe_dir, 'output_phe_' + phe,
                        'HCP_' + phe + '_k_' + str(k) + '_rng_num_' + str(rng))
                    haufe_mat = os.path.join(folder, 'acc.mat')
                    mat = sio.loadmat(haufe_mat)
                    pred_haufe_dict[phe][rng - 1, :] = np.squeeze(
                        mat['optimal_y_p'])
        folder = os.path.join(all_train_dir, 'output_phe_' + phe,
                              'HCP_' + phe + '_k_1019_rng_num_1')
        all_pred_mat = os.path.join(folder, 'acc.mat')
        mat = sio.loadmat(all_pred_mat)
        pred_all_dict[phe] = np.squeeze(mat['optimal_y_p'])

    npy = os.path.join(output_dir, 'HCP_haufe_100.npz')
    np.savez(npy, pred_haufe_dict=pred_haufe_dict)
    npy = os.path.join(output_dir, 'HCP_haufe_all.npz')
    np.savez(npy, pred_all_dict=pred_all_dict)


if __name__ == '__main__':
    get_haufe_pred()
