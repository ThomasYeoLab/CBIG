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
import pandas as pd
import math

from scipy.stats.stats import pearsonr
from sklearn.linear_model import ElasticNetCV
from sklearn.linear_model import ElasticNet

from cbig.CBIG_mics import mics_z_norm, cod_znormed, split_tra_val_with_y
from cbig.config import config


def elasticnet(x_test_k, x_test_remain, y_k, args):
    '''perform elasticnet

    Args:
        x_test_k (ndarray): input data in meta-test set for training
        x_test_remain (ndarray): input data in meta-test set for testing
        y_k (ndarray) : output data in meta-test set for training
        args (argparse.ArgumentParser) : args that could be used by
          other function

    Returns:
        Tuple: prediction for testing and training data in meta-test set

    '''

    parameters = {
        'alpha': [0.00001, 0.0001, 0.001, 0.004, 0.007, 0.01, 0.04, 0.07, 0.1],
        'l1_ratio': [.1, .5, .7, .9, .95, .99, 1]
    }

    regr = ElasticNetCV(l1_ratio=parameters['l1_ratio'],
                        alphas=parameters['alpha'],
                        cv=5,
                        random_state=args.seed,
                        tol=1e-3)
    regr.fit(x_test_k, y_k)

    regr_test = ElasticNet(l1_ratio=regr.l1_ratio_,
                           alpha=regr.alpha_,
                           random_state=args.seed,
                           tol=1e-3)
    regr_test.fit(x_test_k, y_k)

    return regr_test.predict(x_test_remain), regr_test.predict(x_test_k)


def elasticnet_wrapper(x_test, y_test, split_ind, args):
    '''wrapper for elasticnet

    Args:
        x_test (ndarray): input data in meta-test set
        y_test (ndarray): output data in meta-test set
        split_ind (ndarray): training and test split
        args (argparse.ArgumentParser) : args that could be used by
          other function

    Returns:
        res_cor (float) : correlation between prediction and ground truth
        res_cod (float) : COD between prediction and ground truth
        y_pred (ndarray) : prediction for all subjects in meta-test set

    '''

    # get split
    split = np.squeeze(split_ind)
    split_k = split == 0
    split_tes = split == 1
    split_tra, split_val = split_tra_val_with_y(split, y_test)

    assert np.array_equal(split_k, split_tra + split_val)
    _, x_test, _, _ = mics_z_norm(x_test[split_k], x_test)
    x_test_k = x_test[split_k]
    x_test_remain = x_test[split_tes]

    # z norm based on y_train
    _, y_test, _, t_sigma = mics_z_norm(y_test[split_k], y_test)
    y_test_k = y_test[split_k]
    y_test_remain = y_test[split_tes]

    # exclude nan value from y_test_remain
    real_index = ~np.isnan(y_test_remain)
    y_test_remain = y_test_remain[real_index]
    x_test_remain = x_test_remain[real_index]

    y_pred_stack, y_pred_stack_k = elasticnet(x_test_k, x_test_remain,
                                              y_test_k, args)

    res_cor = pearsonr(y_test_remain, y_pred_stack)[0]
    res_cod = cod_znormed(y_test_remain, y_pred_stack)
    y_pred = np.concatenate((y_pred_stack_k, y_pred_stack), axis=0)
    if math.isnan(res_cor):
        raise SystemExit('nan in test correlation')

    return res_cor, res_cod, y_pred


def main(args):
    '''main function for elasticnet

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
    n_rng = args.rng
    dataset = args.dataset

    # load data
    tes_phe = pd.read_csv(os.path.join(args.phe_dir))
    tes_sub_list = pd.read_table(os.path.join(args.sub_dir), header=None)
    tes_mask = tes_phe['eid'].isin(tes_sub_list[0].values.tolist())
    tes_phe_name = tes_phe.columns[args.start_idx:args.end_idx]
    tes_phe = tes_phe.loc[tes_mask, tes_phe_name]
    y_test = tes_phe.values.astype(float)
    x_test = pd.read_csv(os.path.join(args.data_dir))
    x_tes_mask = x_test['eid'].isin(tes_sub_list[0].values.tolist())
    x_test_name = x_test.columns[1:]
    x_test = x_test.loc[x_tes_mask, x_test_name]
    x_test.fillna(x_test.mean(), inplace=True)
    x_test = x_test.values.astype(float)

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
    pred = np.zeros((n_rng, len(ks), len(tes_phe_name), x_test.shape[0]))

    for i in range(n_rng):
        for ik, k in enumerate(ks):
            for ib, phe in enumerate(tes_phe_name):
                split_ind = split_ind_dict.get(phe)[i, ik, :]
                meta_cor[i, ik,
                         ib], meta_cod[i, ik,
                                       ib], tmp_pred = elasticnet_wrapper(
                                           x_test, y_test[:, ib], split_ind,
                                           args)

                pred[i, ik, ib, :tmp_pred.shape[0]] = tmp_pred

    np.savez(meta_cor_npz, meta_cor=meta_cor, meta_cod=meta_cod, pred=pred)

    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
       argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--data_dir', type=str, default=None)
    parser.add_argument('--phe_dir', type=str, default=None)
    parser.add_argument('--sub_dir', type=str, default=None)
    parser.add_argument('--out_dir', '-o', type=str, default=None)
    parser.add_argument('--out_subdir', '-osub', type=str, default=None)
    parser.add_argument('--inter_dir', type=str, default=None)
    parser.add_argument('--in_dir', type=str, default=None)
    parser.add_argument('--dataset', type=str, default=config.DATASET)
    parser.add_argument('--log_stem', type=str, default='elasticnet')
    parser.add_argument('--seed', type=int, default=config.SEED)
    parser.add_argument('--split', type=str, default='test')
    parser.add_argument('--rng', type=int, default=None)
    parser.add_argument("--ks", type=int, nargs='+', default=None)

    fine_tune_parser = parser.add_mutually_exclusive_group(required=False)
    fine_tune_parser.add_argument('--across-dataset',
                                  dest='across_dataset',
                                  action='store_true')
    fine_tune_parser.add_argument('--not-across-dataset',
                                  dest='across_dataset',
                                  action='store_false')
    parser.set_defaults(across_dataset=False)

    parser.add_argument("--start_idx", type=int, default=None)
    parser.add_argument("--end_idx", type=int, default=None)

    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
