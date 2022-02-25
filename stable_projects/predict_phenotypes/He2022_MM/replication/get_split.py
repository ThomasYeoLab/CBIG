#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import argparse
import numpy as np
import scipy.io as sio
from shutil import copyfile


def get_phes(folder, dataset):
    '''get the phenotypes of test meta-set

    Args:
        folder (str): folder path that stores data input files
        dataset (str): state the test meta-set dataset

    Returns:
        List: list of test meta-set phenotypes
    '''

    if dataset == 'ukbb':
        txt = os.path.join(folder, 'ukbb_test_final_phe_list.txt')
    elif dataset == 'HCP':
        txt = os.path.join(folder, 'HCP_diff_roi_final_phe_list.txt')
    else:
        raise NameError('wrong dataset name')
    phes = np.genfromtxt(txt, dtype=str)
    return phes


def get_split(dataset='ukbb'):
    '''save the split files from KRR classical and meta-matching

    Args:
        dataset (str): state the test meta-set dataset

    Returns:
        None
    '''

    base_dir = os.path.join(
        os.environ['CBIG_CODE_DIR'],
        'stable_projects/predict_phenotypes/He2022_MM/replication')

    if dataset == 'ukbb':
        n_subjects = 10000
        classical_dir = os.path.join(base_dir, 'output_KRR_classical_ukbb')
        mm_dir = os.path.join(base_dir, 'output_KRR_mm')
    elif dataset == 'HCP':
        n_subjects = 1019
        classical_dir = os.path.join(base_dir, 'output_KRR_classical_HCP')
    else:
        raise NameError('wrong dataset name')

    output_dir = os.path.join(base_dir, 'output_intermediate')
    os.makedirs(output_dir, exist_ok=True)

    input_dir = os.path.join(os.environ['CBIG_REPDATA_DIR'],
                             'stable_projects/predict_phenotypes/He2022_MM')
    phes_tes = get_phes(input_dir, dataset)
    rngs = 100
    ks = [10, 20, 50, 100, 200]
    split_ind_dict = {}

    for _, phe in enumerate(phes_tes):
        split_ind_dict[phe] = np.zeros((rngs, len(ks), n_subjects))

    # get split file for k shot
    for _, phe in enumerate(phes_tes):
        for ik, k in enumerate(ks):
            for rng in range(1, rngs + 1):
                folder = os.path.join(
                    classical_dir, 'output_phe_' + phe, dataset + '_' + phe +
                    '_k_' + str(k) + '_rng_num_' + str(rng))
                split_mat = os.path.join(folder,
                                         dataset + '_subject_split.mat')
                # print(split_mat)
                mat = sio.loadmat(split_mat)
                split_ind_dict[phe][rng - 1, ik, :] = np.squeeze(
                    mat['sub_fold']['fold_index'][0][0])
        print(phe)

    for _, phe in enumerate(phes_tes):
        split_ind_dict[phe] = split_ind_dict[phe].astype(bool)

    npy = os.path.join(output_dir,
                       dataset + '_split_ind_' + str(rngs) + '.npz')
    np.savez(npy, split_ind_dict=split_ind_dict)

    # get split file from KRR base
    if dataset == 'ukbb':
        rngs = 1
        for rng in range(1, rngs + 1):
            folder = os.path.join(mm_dir, dataset + '_rng_num_' + str(rng))
            split_mat = os.path.join(folder, dataset + '_subject_split.mat')
            target_mat = os.path.join(output_dir,
                                      'split_rng' + str(rng) + '.mat')
            copyfile(split_mat, target_mat)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset", type=str, default="ukbb")
    args = parser.parse_args()
    get_split(args.dataset)
