#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import random
import argparse
import numpy as np
from scipy.stats.stats import pearsonr
import pandas as pd

from sklearn.model_selection import KFold, GridSearchCV
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import normalize

from cbig.CBIG_mics import mics_z_norm, cod_znormed, split_tra_val_with_y
from cbig.config import config


def demean_norm(val):
    '''de-mean and normalize data

    Args:
        val (ndarray): value to be de-meaned and normalized

    Returns:
        ndarray: de-meaned and normalized data
    '''
    mu = np.nanmean(val, axis=1, keepdims=True)
    val = val - mu
    return normalize(val, axis=1, norm='l2')


def stacking(y_pred_k, y_pred_test, y_k, args):
    '''perform KRR for meta-matching stacking

    Args:
        y_pred_k (ndarray): input data in meta-test set for training
        y_pred_test (ndarray): input data in meta-test set for testing
        y_k (ndarray) : output data in meta-test set for training
        args (argparse.ArgumentParser) : args that could be used by
          other function

    Returns:
        Tuple: prediction for testing and training data in meta-test set

    '''

    parameters = {
        'alpha': [5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    }

    krr = KernelRidge()
    cv = KFold(n_splits=5, shuffle=True, random_state=args.seed)
    clf = GridSearchCV(krr, parameters, cv=cv)
    clf.fit(y_pred_k, y_k)

    return clf.predict(y_pred_test), clf.predict(y_pred_k)


def mm_stacking_wrapper(y_test, y_pred, split_ind, args, k):
    '''wrapper for meta-matching stacking

    Args:
        y_test (ndarray): output data in meta-test set
        y_pred (ndarray): prediction in meta-test set from
          pretrained model as input for KRR
        split_ind (ndarray): training and test split
        args (argparse.ArgumentParser) : args that could be used by
          other function
        k (int) : number of k-shot subjects

    Returns:
        res_cor (float) : correlation between prediction and ground truth
        res_cod (float) : COD between prediction and ground truth
        y_pred_new (ndarray) : prediction for all subjects in meta-test set

    '''

    # get split from Classical KRR
    split = np.squeeze(split_ind)
    split_k = split == 0
    split_tes = split == 1
    split_tra, split_val = split_tra_val_with_y(split, y_test)

    assert np.array_equal(split_k, split_tra + split_val)
    y_pred_k = y_pred[split_k, :]
    y_pred_remain = y_pred[split_tes, :]

    # z norm based on y_train
    _, y_test, _, t_sigma = mics_z_norm(y_test[split_k], y_test)
    y_test_k = y_test[split_k]
    y_test_remain = y_test[split_tes]

    # perform stacking
    k_limit = args.k_limit
    dim = min(k, k_limit)
    scores = []
    if args.metric == 'cod':
        for i in range(y_pred.shape[1]):
            tmp = cod_znormed(y_test_k, y_pred_k[:, i])
            y_pred_k_tmp = -y_pred_k[:, i]
            tmp1 = cod_znormed(y_test_k, y_pred_k_tmp)
            if tmp1 > tmp:
                tmp = tmp1
            scores.append(tmp)
    else:
        for i in range(y_pred.shape[1]):
            tmp = pearsonr(y_test_k, y_pred_k[:, i])[0]
            y_pred_k_tmp = -y_pred_k[:, i]
            tmp1 = pearsonr(y_test_k, y_pred_k_tmp)[0]
            if tmp1 > tmp:
                tmp = tmp1
            scores.append(tmp)

    # normalize dnn outputs across train subjects
    _, y_pred, _, _ = mics_z_norm(y_pred[split_k], y_pred)
    y_pred_k = y_pred[split_k]
    y_pred_remain = y_pred[split_tes]

    # exclude nan value from y_test_remain
    real_index = ~np.isnan(y_test_remain)
    y_test_remain = y_test_remain[real_index]
    y_pred_remain = y_pred_remain[real_index, :]

    index = np.array(scores).argsort()[-dim:]
    y_pred_stack, y_pred_stack_k = stacking(y_pred_k[:, index],
                                            y_pred_remain[:, index], y_test_k,
                                            args)

    res_cor = pearsonr(y_test_remain, y_pred_stack)[0]
    res_cod = cod_znormed(y_test_remain, y_pred_stack)
    y_pred_new = np.concatenate((y_pred_stack_k, y_pred_stack), axis=0)

    return res_cor, res_cod, y_pred_new


def main(args):
    '''main function for meta-matching stacking

    Args:
        args: args from command line

    Returns:
        None
    '''

    print('\nCBIG meta-matching (DNN stacking) with argument: ' + str(args))

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)

    log_str = args.log_stem + '_result_' + args.split
    meta_cor_npz = os.path.join(args.out_dir + args.out_subdir,
                                log_str + '.npz')
    os.makedirs(os.path.join(args.out_dir + args.out_subdir), exist_ok=True)
    ks = args.ks
    n_rng = args.n_rng
    dataset = args.dataset

    # load data
    tes_phe = pd.read_csv(os.path.join(args.phe_dir))
    tes_sub_list = pd.read_table(os.path.join(args.sub_dir), header=None)
    tes_mask = tes_phe['eid'].isin(tes_sub_list[0].values.tolist())
    tes_phe_name = tes_phe.columns[args.start_idx:args.end_idx]
    tes_phe = tes_phe.loc[tes_mask, tes_phe_name]
    y_test = tes_phe.values.astype(float)

    # load base prediction result
    npz = os.path.join(args.data_dir, 'dnn_test_base.npz')
    npz = np.load(npz)
    tes_res_record = npz['tes_res_record']

    y_pred = tes_res_record[0, 0, :, :]

    # load data split
    if args.across_dataset:
        npz = np.load(os.path.join(
            args.inter_dir, dataset + '_split_ind_' + str(n_rng) + '.npz'),
                      allow_pickle=True)
    else:
        npz = np.load(os.path.join(args.inter_dir,
                                   'ukbb_split_ind_' + str(n_rng) + '.npz'),
                      allow_pickle=True)
    split_ind_dict = npz['split_ind_dict'].item()

    meta_cor = np.zeros((n_rng, len(ks), len(tes_phe_name)))
    meta_cod = np.zeros((n_rng, len(ks), len(tes_phe_name)))
    pred = np.zeros((n_rng, len(ks), len(tes_phe_name), y_test.shape[0]))

    for i in range(n_rng):
        for ik, k in enumerate(ks):
            for ib, phe in enumerate(tes_phe_name):
                split_ind = split_ind_dict.get(phe)[i, ik, :]
                meta_cor[i, ik,
                         ib], meta_cod[i, ik,
                                       ib], tmp_pred = mm_stacking_wrapper(
                                           y_test[:, ib], y_pred, split_ind,
                                           args, k)

                pred[i, ik, ib, :tmp_pred.shape[0]] = tmp_pred

    np.savez(meta_cor_npz, meta_cor=meta_cor, meta_cod=meta_cod, pred=pred)

    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by
          other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--phe_dir', type=str, default=None)
    parser.add_argument('--dataset', type=str, default=config.DATASET)
    parser.add_argument('--sub_dir', type=str, default=None)
    parser.add_argument('--out_dir', '-o', type=str, default=None)
    parser.add_argument('--out_subdir', '-osub', type=str, default=None)
    parser.add_argument('--inter_dir', type=str, default=None)
    parser.add_argument('--in_dir', type=str, default=None)
    parser.add_argument('--log_stem', type=str, default='meta_stacking')
    parser.add_argument('--seed', type=int, default=config.SEED)
    parser.add_argument('--split', type=str, default='test')
    parser.add_argument('--n_rng', type=int, default=None)
    parser.add_argument('--k_limit', type=int, default=config.K_LIMIT)
    parser.add_argument("--ks", type=int, nargs='+', default=None)

    parser.add_argument('--metric', type=str, default=None)
    parser.add_argument("--start_idx", type=int, default=None)
    parser.add_argument("--end_idx", type=int, default=None)
    parser.add_argument("--data_dir", type=str, default=None)

    across_dataset_parser = parser.add_mutually_exclusive_group(required=False)
    across_dataset_parser.add_argument('--across-dataset',
                                       dest='across_dataset',
                                       action='store_true')
    across_dataset_parser.add_argument('--not-across-dataset',
                                       dest='across_dataset',
                                       action='store_false')
    parser.set_defaults(across_dataset=False)

    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
