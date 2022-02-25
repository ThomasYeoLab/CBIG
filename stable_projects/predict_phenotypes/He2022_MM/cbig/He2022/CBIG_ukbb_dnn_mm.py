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
from scipy.stats.stats import pearsonr

from config import config
from CBIG_mics import mics_z_norm, cod_znormed


def mm_dnn(y_test, y_pred, split_ind):
    '''meta-matching with DNN

    Args:
        y_test (ndarray): original test data
        y_pred (ndarray): predicted for test subjects with base model trained
            on training meta-set
        split_ind (ndarray): array that indicate K subjects

    Returns:
        Tuple: result (correlation and COD) of meta-matching with DNN, best
            matched phenotypes, prediction of meta-matching on k shot subjects
    '''
    # get split from krr pure baseline
    split = np.squeeze(split_ind)
    split_k = split == 0
    split_tes = split == 1

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

    best_score = float("-inf")
    best_phe = -1
    best_sign = 1
    for i in range(y_pred.shape[1]):
        sign = 1
        tmp = cod_znormed(y_test_k, y_pred_k[:, i])
        y_pred_k_tmp = -y_pred_k[:, i]
        tmp1 = cod_znormed(y_test_k, y_pred_k_tmp)
        if tmp1 > tmp:
            tmp = tmp1
            sign = -1
        if tmp >= best_score:
            best_score = tmp
            best_phe = i
            best_sign = sign
    y_mm_pred = best_sign * y_pred_remain[:, best_phe]
    y_mm_pred_k = best_sign * y_pred_k[:, best_phe]
    res_cor = pearsonr(y_test_remain, y_mm_pred)[0]
    res_cod = cod_znormed(y_test_remain, y_mm_pred)

    return res_cor, res_cod, best_phe, y_mm_pred_k


def main(args):
    '''main function for meta-matching (DNN)

    Args:
        args: args from command line

    Returns:
        None
    '''

    print('\nCBIG meta-matching (DNN) with argument: ' + str(args))

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)

    n_rng = 100
    tuned_by = 'cod'
    ks = [10, 20, 50, 100, 200]
    log_str = args.log_stem + '_result_' + args.split
    if args.across_dataset:
        log_str += '_across_dataset'
    meta_cor_npz = os.path.join(args.out_dir, log_str + '.npz')
    meta_temp_npz = os.path.join(args.out_dir, 'tmp/' + log_str + '_rng.npz')
    os.makedirs(os.path.join(args.out_dir, 'tmp'), exist_ok=True)

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
    # x_test = npz['x_test']
    y_test = npz['y_test_different_set_phe']
    tes_phe = npz['tes_phe']
    tra_phe = npz['tra_phe']

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

    # perform meta matching with DNN and transfer learning
    start_time = time.time()
    meta_cor = np.zeros((n_rng, len(ks), len(tes_phe)))
    meta_cod = np.zeros((n_rng, len(ks), len(tes_phe)))
    meta_phe = np.zeros((n_rng, len(ks), len(tes_phe)))
    if args.haufe_save:
        meta_pred_k_100 = np.zeros((n_rng, len(tes_phe), 100))
    for i in range(n_rng):
        for ik, k in enumerate(ks):
            for ib, phe in enumerate(tes_phe):
                split_ind = split_ind_dict.get(phe)[i, ik, :]
                meta_cor[i, ik, ib], meta_cod[
                    i, ik, ib], meta_phe[i, ik, ib], tmp_pred = mm_dnn(
                        y_test[:, ib], y_pred, split_ind)
                if args.haufe_save and k == 100:
                    meta_pred_k_100[i, ib, :] = tmp_pred
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
            meta_temp_npz,
            meta_cor=meta_cor,
            meta_cod=meta_cod,
            current_rng=i,
            meta_phe=meta_phe,
            tes_phe=tes_phe,
            tra_phe=tra_phe)

    np.savez(
        meta_cor_npz,
        meta_cor=meta_cor,
        meta_cod=meta_cod,
        meta_phe=meta_phe,
        tes_phe=tes_phe,
        tra_phe=tra_phe)
    if args.haufe_save:
        haufe_npz = os.path.join(args.large_data_dir,
                                 'haufe_y_pred_100_basic_mm.npz')
        np.savez(haufe_npz, meta_pred_k_100=meta_pred_k_100)

    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--out_dir', '-o', type=str, default=config.OUT_DIR)
    parser.add_argument('--in_dir', type=str, default=config.IN_DIR)
    parser.add_argument('--inter_dir', type=str, default=config.INTER_DIR)
    parser.add_argument(
        '--large_data_dir', type=str, default=config.LARGE_DATA_DIR)
    parser.add_argument('--log_stem', type=str, default='meta')
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--split', type=str, default='test')

    across_dataset_parser = parser.add_mutually_exclusive_group(required=False)
    across_dataset_parser.add_argument(
        '--across-dataset', dest='across_dataset', action='store_true')
    across_dataset_parser.add_argument(
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
