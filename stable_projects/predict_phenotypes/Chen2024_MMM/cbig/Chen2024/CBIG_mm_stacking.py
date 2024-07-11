#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""
import os
import time
import random
import argparse
import numpy as np
from scipy.stats.stats import pearsonr
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import KFold, GridSearchCV
from config import config
from CBIG_misc import read_phes, get_phe_num, misc_z_norm, cod_znormed, split_tra_val_with_y, demean_normalize


def KRR_stacking(y_pred_k, y_pred_test, y_k, args):
    '''perform stacking using KRR model

    Args:
        y_pred_k (ndarray): prediction for K subjects in meta-test dataset
        y_pred_test (ndarray): prediction for remaining meta-test subjects
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
    cv = KFold(n_splits=config.N_CV, shuffle=True, random_state=args.seed)
    clf = GridSearchCV(krr, parameters, cv=cv)
    clf.fit(demean_normalize(y_pred_k), y_k)

    return clf.predict(demean_normalize(y_pred_test)), clf.predict(demean_normalize(y_pred_k))


def mm_stacking(y_test, y_pred, split_ind, args):
    '''meta-matching with stacking

    Args:
        y_pred_k (ndarray): prediction for K subjects in meta-test dataset
        y_pred_test (ndarray): prediction for remaining meta-test subjects
        split_ind (ndarray): array that indicate K subjects
        args: args from command line

    Returns:
        Tuple: result (correlation) of meta-matching with stacking,
            result (COD) of meta-matching with stacking, prediction
            of stacking on k shot subjects, split index of k shot subjects.
    '''

    # get split from Classical KRR
    split = np.squeeze(split_ind)
    split_k = split == False
    split_tes = split == True
    split_tra, split_val = split_tra_val_with_y(split, y_test)

    assert np.array_equal(split_k, split_tra + split_val)
    y_pred_k = y_pred[split_k, :]
    y_pred_remain = y_pred[split_tes, :]

    # z norm based on K-shot (training set)
    _, y_test, _, _ = misc_z_norm(y_test[split_k], y_test)

    y_test_k = y_test[split_k]
    y_test_remain = y_test[split_tes]

    # exclude nan value from y_test_remain
    real_index = ~np.isnan(y_test_remain)
    y_test_remain = y_test_remain[real_index]
    y_pred_remain = y_pred_remain[real_index, :]

    # perform stacking using KRR (we use all source phenotypic prediction rather than choosing top-K)
    y_pred_stack, y_pred_stack_k = KRR_stacking(y_pred_k, y_pred_remain, y_test_k, args)

    res_cor = pearsonr(y_test_remain, y_pred_stack)[0]
    res_cod = cod_znormed(y_test_remain, y_pred_stack)

    return res_cor, res_cod, y_pred_stack_k, split_k


def main(args):
    '''main function for CBIG meta-matching with stacking

    Args:
        args: args from command line

    Returns:
        None
    '''
    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    ks = config.KS
    n_rng = config.N_RNG


    if args.exp_dataset:
        args.in_dir = config.IN_DIR_UT if args.unit_test else config.IN_DIR_EXP
        args.out_dir = config.OUT_DIR_UT if args.unit_test else config.OUT_DIR_EXP
        args.inter_dir = config.INTER_DIR_UT if args.unit_test else config.INTER_DIR_EXP
        args.model_dir = config.MODEL_DIR_UT if args.unit_test else config.MODEL_DIR_EXP
        args.tar_dataset = config.DATASET_NAME_EXP['test'][0]
        args.dataset_name = config.DATASET_NAME_EXP
        n_rng = config.N_RNG_EXP

    print('\nCBIG meta-matching (stacking) with argument: ' + str(args))

    tar_dataset = args.tar_dataset

    log_str = args.log_stem + '_2' + tar_dataset + '_result'
    meta_cor_npz = os.path.join(args.out_dir, log_str + '.npz')
    os.makedirs(os.path.join(args.out_dir), exist_ok=True)

    # load base prediction result
    dnn_pred = np.load(os.path.join(args.inter_dir, 'dnn_prediction.npz'))
    y_pred_dnn_xlarge = dnn_pred[tar_dataset]
    n_test_subj = y_pred_dnn_xlarge.shape[0]

    if "MM_stacking" in args.log_stem:
        y_pred = y_pred_dnn_xlarge
    elif "dataset_stacking" in args.log_stem:
        rr_pred = np.load(os.path.join(args.inter_dir, 'rr_prediction.npz'), allow_pickle=True)
        y_pred_rr_xlarge = rr_pred[tar_dataset]
        y_pred_rr_1layer = {}
        for src_dataset in args.dataset_name['large'] + args.dataset_name['medium']:
            n_phe = get_phe_num(args.in_dir, src_dataset)
            y_pred_rr_1layer[src_dataset] = np.zeros((n_test_subj, n_phe))
            for i in range(n_phe):
                npz = np.load(os.path.join(args.inter_dir, src_dataset + '_rr_base_' + str(i) + '.npz'))
                y_pred_rr_1layer[src_dataset][:, i] = npz[tar_dataset].reshape((n_test_subj,))
        y_pred = np.concatenate([y_pred_dnn_xlarge] + [y_pred_rr_xlarge] + list(y_pred_rr_1layer.values()), axis=1)
    elif 'multilayer_stacking' in args.log_stem:
        rr_pred = np.load(os.path.join(args.inter_dir, 'rr_prediction.npz'), allow_pickle=True)
        y_pred_rr_xlarge = rr_pred[tar_dataset]
        y_pred_rr_1layer = {}
        y_pred_rr_2layer = {}
        for src_dataset in args.dataset_name['large'] + args.dataset_name['medium']:
            n_phe = get_phe_num(args.in_dir, src_dataset)
            y_pred_rr_1layer[src_dataset] = np.zeros((n_test_subj, n_phe))
            y_pred_rr_2layer[src_dataset] = np.zeros((n_test_subj, n_phe))
            for i in range(n_phe):
                base_npz = np.load(os.path.join(args.inter_dir, src_dataset + '_rr_base_' + str(i) + '.npz'))
                y_pred_rr_1layer[src_dataset][:, i] = base_npz[tar_dataset].reshape((n_test_subj,))

                multilayer_npz = np.load(os.path.join(
                    args.inter_dir, src_dataset + '_rr_multilayer_' + str(i) + '.npz'))
                y_pred_rr_2layer[src_dataset][:, i] = multilayer_npz[tar_dataset].reshape((n_test_subj,))
        y_pred = np.concatenate([y_pred_dnn_xlarge] + [y_pred_rr_xlarge] +
                                 list(y_pred_rr_1layer.values()) +
                                 list(y_pred_rr_2layer.values()), axis=1)
    else:
        raise NameError("Wrong method name (try 'MM_stacking', 'dataset_stacking', 'multilayer_stacking')")

    # load original data
    print(y_pred.shape)
    dnn_input = np.load(os.path.join(args.in_dir, tar_dataset, tar_dataset + '_dnn_input.npz'), allow_pickle=True)
    y_test = dnn_input['y_raw']
    y_test = np.array(y_test, dtype=np.float64)

    tes_phe = read_phes(args.in_dir, tar_dataset)

    if args.haufe_save:
        x_test = dnn_input['x_raw']

    # load data split
    split = np.load(
        os.path.join(args.inter_dir, tar_dataset + '_split_ind_' + str(n_rng) + '.npz'),
        allow_pickle=True)
    split_ind_dict = split['split_ind_dict'].item()

    # perform meta matching with stacking
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
                    y_test[:, ib], y_pred, split_ind, args)
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
        meta_cor_npz, meta_cor=meta_cor, meta_cod=meta_cod, tes_phe=tes_phe)

    if args.haufe_save:
        haufe_npz = os.path.join(args.large_data_dir,
                                 'haufe_y_pred_100_' + args.log_stem + '_2' + tar_dataset + '.npz')
        np.savez(haufe_npz, meta_pred_k_100=meta_pred_k_100, meta_x_k_100=meta_x_k_100)

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
    parser.add_argument('--log_stem', type=str, default='multilayer_stacking')
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--tar_dataset', type=str, default='HCP')
    parser.add_argument('--dataset_name', type=dict, default=config.DATASET_NAME)

    restricted_alpha_parser = parser.add_mutually_exclusive_group(required=False)
    restricted_alpha_parser.add_argument('--restricted-alpha', dest='restricted_alpha', action='store_true')
    restricted_alpha_parser.add_argument('--not-restricted-alpha', dest='restricted_alpha', action='store_false')
    parser.set_defaults(restricted_alpha=True)

    haufe_save_parser = parser.add_mutually_exclusive_group(required=False)
    haufe_save_parser.add_argument('--haufe-save', dest='haufe_save', action='store_true')
    haufe_save_parser.add_argument('--not-haufe-save', dest='haufe_save', action='store_false')
    parser.set_defaults(haufe_save=False)

    exp_dataset_parser = parser.add_mutually_exclusive_group(required=False)
    exp_dataset_parser.add_argument('--exp-dataset', dest='exp_dataset', action='store_true')
    exp_dataset_parser.add_argument('--not-exp-dataset', dest='exp_dataset', action='store_false')
    parser.set_defaults(exp_dataset=False)

    unit_tests_parser = parser.add_mutually_exclusive_group(required=False)
    unit_tests_parser.add_argument('--unit-test', dest='unit_test', action='store_true')
    unit_tests_parser.add_argument('--not-unit-test', dest='unit_test', action='store_false')
    parser.set_defaults(unit_test=False)
    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
