#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import argparse
import numpy as np
import os
import pandas as pd
import random
from cbig.config import config


def rearrange_hist_center(y_train):
    '''function to get centers for histogram of y for training subjects

    Args:
       y_train (ndarray): y (output) data of training subjects in
                         meta-test data

    Returns:
       hist_centers (ndarray) : centers for histogram of y

    '''

    tmp_hist_edges = np.unique(y_train)
    if len(tmp_hist_edges) == 1:
        hist_centers = tmp_hist_edges
    else:
        interval = (tmp_hist_edges[1:] - tmp_hist_edges[:-1]) / 2
        tmp_hist_centers = list(tmp_hist_edges[:-1] + interval)
        hist_centers = np.hstack([tmp_hist_centers[0] - interval[0] * 2] +
                                 tmp_hist_centers +
                                 [tmp_hist_centers[-1] + interval[-1] * 2])

    return hist_centers


def split_sub(seed, num_real, sub_list, subject_list_real, subject_list_nan,
              num_total, num_test):
    '''function to get subject index for training and testing

    Args:
        seed (int): random seed
        num_real (ndarray): number of y with real value
        sub_list (ndarray): subject list
        subject_list_real (ndarray) : subject list with real value of y
        subject_list_nan (ndarray) : subject list with NaN value of y
        num_total (int) : number of total subjects
        num_test (int) : number of test subjects

    Returns:
        fold_index (ndarray) : subject index for training and testing

    '''

    rand_idx = np.random.RandomState(seed=seed).permutation(num_real)
    subject_list_real_tmp = subject_list_real[rand_idx]

    # generate list
    subject_list_rand = np.hstack([subject_list_real_tmp, subject_list_nan])
    subject_list_test = np.array(sorted(subject_list_rand[-num_test:]))
    fold_index = np.zeros((num_total))
    fold_index[np.where(np.isin(sub_list, subject_list_test))[0]] = 1

    return fold_index


def generate_split_kshot_mm(args):
    '''main function for generating k-shot split in meta-test set
        for 100 repetitions and for all phenotypes

    Args:
        args: args from command line

    Returns:
        None
    '''

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    num_inner_folds = args.num_inner_folds
    thr_y_znorm = 1.0 - 1.0 / num_inner_folds

    ks = args.ks
    start_idx = args.start_idx
    end_idx = args.end_idx
    n_rng = args.n_rng
    phe = pd.read_csv(args.phe_dir)
    sub_list = pd.read_table(os.path.join(args.sub_dir), header=None)
    mask = phe['eid'].isin(sub_list[0].values.tolist())
    phe_name = phe.columns[start_idx:end_idx]
    phe = phe.loc[mask, phe_name]
    y_test_different_set_phe = phe.values.astype(float)

    # Train/Validation/Test split
    num_total = len(sub_list)
    sub_list = np.hstack(sub_list.values.tolist())

    split_ind_dict = {}
    for _, phe in enumerate(phe_name):
        split_ind_dict[phe] = np.zeros((n_rng, len(ks), num_total))

    for ib, phe in enumerate(phe_name):
        y = y_test_different_set_phe[:, ib]

        for ik, k in enumerate(ks):
            num_train = k
            num_test = num_total - num_train
            # handle missing value
            ind_real = ~np.isnan(y)
            num_real = sum(ind_real)
            subject_list_real = sub_list[ind_real]
            subject_list_nan = sub_list[~ind_real]
            rng_offset = 0

            for i in range(n_rng):
                seed = i + rng_offset
                fold_index = split_sub(seed, num_real, sub_list,
                                       subject_list_real, subject_list_nan,
                                       num_total, num_test)

                y_train = y[fold_index == 0]
                hist_centers = rearrange_hist_center(y_train)

                while (np.unique(y_train).shape[0] == 1) or \
                    (max(np.histogram(y_train, hist_centers)[0])
                     / y_train.shape[0] >= thr_y_znorm):
                    print('rng ' + str(i + rng_offset) + ' for ' + phe +
                          ' has ' + str(thr_y_znorm) + ' same value')
                    rng_offset = rng_offset + 1

                    seed = i + rng_offset
                    fold_index = split_sub(seed, num_real, sub_list,
                                           subject_list_real, subject_list_nan,
                                           num_total, num_test)
                    y_train = y[fold_index == 0]
                    hist_centers = rearrange_hist_center(y_train)

                split_ind_dict[phe][i, ik, :] = fold_index

    for _, phe in enumerate(phe_name):
        split_ind_dict[phe] = split_ind_dict[phe].astype(bool)

    if not os.path.exists(args.inter_dir):
        os.makedirs(args.inter_dir)
    npy = os.path.join(args.inter_dir,
                       args.dataset + '_split_ind_' + str(n_rng) + '.npz')
    np.savez(npy, split_ind_dict=split_ind_dict)


def get_args():
    '''function to get args from command line and return the args

    Returns:
       argparse.ArgumentParser: args that could be used by other function
    '''

    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--phe_dir', type=str, default=None)
    parser.add_argument('--sub_dir', type=str, default=None)
    parser.add_argument('--inter_dir', type=str, default=None)
    parser.add_argument('--dataset', type=str, default='ukbb')
    parser.add_argument('--seed', type=int, default=config.SEED)
    parser.add_argument('--n_rng', type=int, default=None)
    parser.add_argument('--start_idx', type=int, default=None)
    parser.add_argument('--end_idx', type=int, default=None)
    parser.add_argument('--num_inner_folds', type=int, default=None)
    parser.add_argument("--ks", type=int, nargs='+', default=None)

    return parser.parse_args()


if __name__ == '__main__':
    generate_split_kshot_mm(get_args())
