#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import time
import random
import argparse
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pickle
import nibabel as nib

from cbig.CBIG_mics import read_datapath
from cbig.CBIG_model_pytorch import crop_center
from cbig.config import config

plt.switch_backend('Agg')


def load_3D_input(sublist):
    '''load volumetric T1 data
       Args:
           sublist (list): list of participants
       Returns:
           data_arr (ndarray): T1 data of all participants in list
    '''

    cutoff_size = (160, 192, 160)

    data_list = []
    for idx in range(len(sublist)):
        nifti_img = nib.load(sublist[idx])
        nii_data = nifti_img.get_fdata()

        nii_data = crop_center(nii_data, cutoff_size)
        nii_data[nii_data < 0] = 0
        x = (nii_data / nii_data.mean()).flatten()
        data_list.append(x)

    data_arr = np.array(data_list)
    return data_arr


def sum_of_mul(A, B):
    '''sum of multiplication of two array over axis=1
    Args:
        A (ndarray): first array for calculation
        B (ndarray): second array for calculation
    Returns:
        ndarray: sum of multiplication calculated
    '''

    return np.einsum('ij,ij->i', A, B)


def covariance_rowwise(A, B):
    '''compute rowwise covariance
    Args:
       A (ndarray): first array for covariance calculation, n_subject x
          n_features
       B (ndarray): second array for covariance calculation, n_subject x 1
    Returns:
       ndarray: covariance calculated
    '''

    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(0, keepdims=True)
    B_mB = B - B.mean(0, keepdims=True)

    N = A.shape[0]
    if B_mB.ndim == 1:
        B_mB = np.expand_dims(B_mB, -1)
    a_nsample = A_mA.shape[1]
    b_nsample = B_mB.shape[1]
    rnt = np.zeros((a_nsample, b_nsample))
    comb = np.array(list(itertools.product(range(a_nsample),
                                           range(b_nsample))))

    n_comb = len(comb)
    chunk = 100000
    if n_comb > chunk:
        start_time = time.time()
        cov = np.empty(n_comb)
        for i in range(chunk, n_comb, chunk):
            cov[i - chunk:i] = sum_of_mul(A_mA[:, comb[i - chunk:i, 0]].T,
                                          B_mB[:, comb[i - chunk:i, 1]].T)

            print(i, time.time() - start_time)
        cov[i:] = sum_of_mul(A_mA[:, comb[i:, 0]].T, B_mB[:, comb[i:, 1]].T)
    else:
        cov = sum_of_mul(A_mA[:, comb[:, 0]].T, B_mB[:, comb[:, 1]].T)
    rnt[comb[:, 0], comb[:, 1]] = cov
    return np.squeeze(rnt) / (N - 1)


def compute_PNF(x, y_pred, ik=None, split_ind_dict=None, phe=None, ib=None):
    '''compute predictive network features (PNF) with various shape
    Args:
        x (ndarray): array of FC for PNF calculation
        y_pred (ndarray): array of predicted y for PNF calculation
        ik (int): index of K-shot used for calculate haufe
        split_ind_dict (dictionary): split file
        phe (str): phenotype name used for calculate haufe
        ib (int): index of phenotype used for calculate haufe

    Returns:
        ndarray: PNF calculated
    '''

    if len(y_pred.shape) < 3:  # for 1 split
        cov = covariance_rowwise(x, y_pred)
        print(cov.shape)
    else:
        cov = np.zeros((100, x.shape[-1]))  # for 100 random splits
        start_time = time.time()

        for i in range(100):
            print(i, time.time() - start_time)

            split_ind = split_ind_dict.get(phe)[i, ik, :]
            split = np.squeeze(split_ind)
            split_k = split == 0

            cov[i, :] = covariance_rowwise(x[split_k, :], y_pred[i, ib, :])
    return cov


def haufe_transform_check(args):
    '''main function for Haufe transform check
    Args:
        args: args from command line
    Returns:
        None
    '''

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    ik = args.ik
    k = args.k
    n_rng = args.rng
    dataset = args.dataset
    phe = args.phe
    ib = args.ib

    # load original data
    tes_sub_list = pd.read_table(os.path.join(args.sub_dir), header=None)
    x_test_path_list = read_datapath(args.data_dir,
                                     list(map(str, tes_sub_list[0])),
                                     args.dataset)
    x_test = load_3D_input(x_test_path_list)

    npz = np.load(os.path.join(args.inter_dir,
                               dataset + '_split_ind_' + str(n_rng) + '.npz'),
                  allow_pickle=True)
    split_ind_dict = npz['split_ind_dict'].item()

    os.makedirs(os.path.join(args.out_dir, 'haufe_check_results'),
                exist_ok=True)
    res_npz_file = os.path.join(
        args.out_dir, 'haufe_check_results',
        'haufe_res_PNF_' + args.stem + '_' + phe + '.npz')
    res_npz = {}

    file_location = 'mm_stacking/meta_stacking_result_test.npz'

    npz = os.path.join(args.out_dir, file_location)
    npz = np.load(npz)
    y_pred = npz['pred'][:, ik, :, :k]

    pnf = compute_PNF(x_test, y_pred, ik, split_ind_dict, phe, ib)
    res_npz['pnf'] = pnf

    with open(res_npz_file, 'wb') as outfile:
        pickle.dump(res_npz, outfile, pickle.HIGHEST_PROTOCOL)

    return


def get_args():
    '''function to get args from command line and return the args
    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--sub_dir', type=str, default=None)
    parser.add_argument('--data_dir', type=str, default=None)
    parser.add_argument('--inter_dir', type=str, default=None)
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--dataset', type=str, default=None)
    parser.add_argument('--stem', type=str, default=None)
    parser.add_argument('--ik', type=int, default=None)
    parser.add_argument('--k', type=int, default=None)
    parser.add_argument('--ib', type=int, default=None)
    parser.add_argument('--phe', type=str, default=None)
    parser.add_argument('--rng', type=int, default=None)
    parser.add_argument('--seed', type=int, default=config.SEED)

    return parser.parse_args()


if __name__ == '__main__':
    haufe_transform_check(get_args())
