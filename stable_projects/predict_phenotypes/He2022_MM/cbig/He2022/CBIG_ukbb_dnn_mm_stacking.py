#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import time
import random
import argparse
import numpy as np
from config import config
from scipy.stats.stats import pearsonr
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import normalize
from CBIG_mics import mics_z_norm, cod_znormed, split_tra_val_with_y


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
    '''perform stacking

    Args:
        y_pred_k (ndarray): predicted for K subjects with base model trained
            on training meta-set
        y_pred_k (ndarray): predicted for remaining test subjects with base
            model trained on training meta-set
        y_k (ndarray): original test data on k subjects
        args: args from command line

    Returns:
        Tuple: predicted value on remaining test subjects and k shot subjects
            with stacking
    '''
    if args.restricted_alpha:
        parameters = {
            'alpha': [5, 10, 15, 20],
        }
    else:
        parameters = {
            'alpha': [
                0.00001, 0.0001, 0.001, 0.004, 0.007, 0.01, 0.04, 0.07, 0.1,
                0.4, 0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 15, 20
            ],
        }
    krr = KernelRidge()
    cv = KFold(n_splits=5, shuffle=True, random_state=args.seed)
    clf = GridSearchCV(krr, parameters, cv=cv)
    clf.fit(demean_norm(y_pred_k), y_k)
    return clf.predict(demean_norm(y_pred_test)), clf.predict(
        demean_norm(y_pred_k))


def mm_stacking(y_test, y_pred, k, split_ind, args):
    '''meta-matching with DNN and stacking

    Args:
        y_test (ndarray): original test data
        y_pred (ndarray): predicted for test subjects with base model trained
            on training meta-set
        k (int): number of k subjects
        split_ind (ndarray): array that indicate K subjects
        args: args from command line

    Returns:
        Tuple: result (correlation) of meta-matching with DNN and stacking,
            result (COD) of meta-matching with DNN and stacking, prediction
            of stacking on k shot subjects, split index of k shot subjects.
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
    _, y_test, _, _ = mics_z_norm(y_test[split_k], y_test)
    y_test_k = y_test[split_k]
    y_test_remain = y_test[split_tes]

    # exclude nan value from y_test_remain
    real_index = ~np.isnan(y_test_remain)
    y_test_remain = y_test_remain[real_index]
    y_pred_remain = y_pred_remain[real_index, :]

    # perform stacking
    k_limit = args.k_limit
    dim = min(k, k_limit)
    scores = []
    for i in range(y_pred.shape[1]):
        tmp = cod_znormed(y_test_k, y_pred_k[:, i])
        y_pred_k_tmp = -y_pred_k[:, i]
        tmp1 = cod_znormed(y_test_k, y_pred_k_tmp)
        if tmp1 > tmp:
            tmp = tmp1
        scores.append(tmp)
    index = np.array(scores).argsort()[-dim:]
    y_pred_stack, y_pred_stack_k = stacking(
        y_pred_k[:, index], y_pred_remain[:, index], y_test_k, args)

    res_cor = pearsonr(y_test_remain, y_pred_stack)[0]
    res_cod = cod_znormed(y_test_remain, y_pred_stack)

    return res_cor, res_cod, y_pred_stack_k, split_k


def main(args):
    '''main function for CBIG meta-matching (DNN stacking)

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

    tuned_by = 'cod'
    log_str = args.log_stem + '_result_' + args.split
    meta_cor_npz = os.path.join(args.out_dir, log_str + '.npz')
    os.makedirs(os.path.join(args.out_dir), exist_ok=True)
    ks = [10, 20, 50, 100, 200]
    n_rng = 100

    # load base prediction result
    if args.across_dataset:
        npz = os.path.join(args.out_dir, 'dnn_across_dataset_base.npz')
    else:
        npz = os.path.join(args.out_dir, 'dnn_base.npz')
    npz = np.load(npz)
    tes_res_record = npz['tes_res_record']
    val_record = npz['val_' + tuned_by + '_record']
    temp = np.mean(val_record[0, :, :], axis=1)
    temp = np.convolve(temp, np.ones(3, dtype=int), 'valid') / 3
    index = np.nanargmax(temp)
    index = index + 1
    print('\nBest validation at index: ', index)
    y_pred = tes_res_record[0, index, :, :]

    # load original data
    if args.across_dataset:
        npz = os.path.join(args.in_dir, 'ukbb_dnn_input_cross_dataset.npz')
    else:
        npz = os.path.join(args.in_dir, 'ukbb_dnn_input_test.npz')
    npz = np.load(npz, allow_pickle=True)
    y_test = npz['y_test_different_set_phe']
    tes_phe = npz['tes_phe']
    tra_phe = npz['tra_phe']
    if args.haufe_save:
        x_test = npz['x_test']

    # load data split
    if args.across_dataset:
        npz = np.load(
            os.path.join(args.inter_dir,
                         'HCP_split_ind_' + str(n_rng) + '.npz'),
            allow_pickle=True)
    else:
        npz = np.load(
            os.path.join(args.inter_dir,
                         'ukbb_split_ind_' + str(n_rng) + '.npz'),
            allow_pickle=True)
    split_ind_dict = npz['split_ind_dict'].item()

    # perform meta matching with FNN and stacking
    start_time = time.time()
    meta_cor = np.zeros((n_rng, len(ks), len(tes_phe)))
    meta_cod = np.zeros((n_rng, len(ks), len(tes_phe)))
    if args.haufe_save:
        meta_pred_k_100 = np.zeros((n_rng, len(tes_phe), 100))
        meta_x_k_100 = np.zeros((n_rng, len(tes_phe), 100, x_test.shape[1]))
    for i in range(n_rng):
        for ik, k in enumerate(ks):
            for ib, phe in enumerate(tes_phe):
                split_ind = split_ind_dict.get(phe)[i, ik, :]
                meta_cor[i, ik, ib], meta_cod[
                    i, ik, ib], tmp_pred, split_k = mm_stacking(
                        y_test[:, ib], y_pred, k, split_ind, args)
                if args.haufe_save and k == 100:
                    meta_pred_k_100[i, ib, :] = tmp_pred
                    meta_x_k_100[i, ib, :, :] = x_test[split_k, :]
        print("rng %d at %ss: cor %.5f, cod %.5f" %
              (i, time.time() - start_time, np.nanmean(meta_cor[:i + 1, :, :]),
               np.nanmean(meta_cod[:i + 1, :, :])))
        mean_cor = np.squeeze(
            np.nanmean(np.nanmean(meta_cor[:i + 1, :, :], axis=2), axis=0))
        mean_cod = np.squeeze(
            np.nanmean(np.nanmean(meta_cod[:i + 1, :, :], axis=2), axis=0))
        print(' '.join('%.6f' % tmp for tmp in mean_cor), ' COD ',
              ' '.join('%.6f' % tmp for tmp in mean_cod))
    np.savez(
        meta_cor_npz,
        meta_cor=meta_cor,
        meta_cod=meta_cod,
        tes_phe=tes_phe,
        tra_phe=tra_phe)

    if args.haufe_save:
        if args.restricted_alpha:
            haufe_npz = os.path.join(
                args.large_data_dir,
                'haufe_y_pred_100_stacking_restricted_alpha.npz')
        else:
            haufe_npz = os.path.join(args.large_data_dir,
                                     'haufe_y_pred_100_stacking.npz')
        np.savez(
            haufe_npz,
            meta_pred_k_100=meta_pred_k_100,
            meta_x_k_100=meta_x_k_100)

    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--out_dir', '-o', type=str, default=config.OUT_DIR)
    parser.add_argument('--inter_dir', type=str, default=config.INTER_DIR)
    parser.add_argument('--in_dir', type=str, default=config.IN_DIR)
    parser.add_argument(
        '--large_data_dir', type=str, default=config.LARGE_DATA_DIR)
    parser.add_argument('--log_stem', type=str, default='meta_stacking')
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--split', type=str, default='test')
    parser.add_argument('--k_limit', type=int, default=33)

    restricted_alpha_parser = parser.add_mutually_exclusive_group(
        required=False)
    restricted_alpha_parser.add_argument(
        '--restricted-alpha', dest='restricted_alpha', action='store_true')
    restricted_alpha_parser.add_argument(
        '--not-restricted-alpha',
        dest='restricted_alpha',
        action='store_false')
    parser.set_defaults(restricted_alpha=True)

    fine_tune_parser = parser.add_mutually_exclusive_group(required=False)
    fine_tune_parser.add_argument(
        '--across-dataset', dest='across_dataset', action='store_true')
    fine_tune_parser.add_argument(
        '--not-across-dataset', dest='across_dataset', action='store_false')
    parser.set_defaults(across_dataset=False)

    haufe_save_parser = parser.add_mutually_exclusive_group(required=False)
    haufe_save_parser.add_argument(
        '--haufe-save', dest='haufe_save', action='store_true')
    haufe_save_parser.add_argument(
        '--not-haufe-save', dest='haufe_save', action='store_false')
    parser.set_defaults(haufe_save=False)

    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
