#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import random
import argparse
import numpy as np
from config import config
from CBIG_rr_large import RR
from CBIG_misc import misc_z_norm, get_phe_num


def main(args):
    '''train Ridge Regression base models on medium source datasets (e.g., GSP/HBN/eNKI) and predict on test datasets

    Args:
        args: args from command line

    Returns:
        None
    '''
    # set all the seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    if args.exp_dataset:
        args.in_dir = config.IN_DIR_UT if args.unit_test else config.IN_DIR_EXP
        args.out_dir = config.OUT_DIR_UT if args.unit_test else config.OUT_DIR_EXP
        args.inter_dir = config.INTER_DIR_UT if args.unit_test else config.INTER_DIR_EXP
        args.model_dir = config.MODEL_DIR_UT if args.unit_test else config.MODEL_DIR_EXP
        args.src_dataset = config.DATASET_NAME_EXP['medium'][0]
        args.dataset_name = config.DATASET_NAME_EXP

    phe_idx = args.phe_idx
    src_dataset = args.src_dataset

    # load input and label data
    if args.multilayer == False:
        # apply direct meta-matching
        npz = np.load(os.path.join(args.in_dir, src_dataset, src_dataset + '_dnn_input.npz'), allow_pickle=True)
        x_train = npz['x_raw']
        x_train[np.isnan(x_train)] = 0
        y_train = npz['y_raw']
        y_train = np.array(y_train, dtype=np.float64)

        x_test = {}  # read FC data in all test datasets
        for tar_dataset in args.dataset_name['test']:
            x_test[tar_dataset] = np.load(os.path.join(args.in_dir, tar_dataset, tar_dataset +
                                                         '_dnn_input.npz'), allow_pickle=True)['x_raw']
            x_test[tar_dataset][np.isnan(x_test[tar_dataset])] = 0
    else:
        # perform multi-layer meta-matching (2-layer)
        # prediction from DNN and RR models trained on extra-large source dataset
        dnn_pred = np.load(os.path.join(args.inter_dir, 'dnn_prediction.npz'))
        x_train_dnn_pred = dnn_pred[src_dataset]
        x_test_dnn_pred = {}
        for name_test in args.dataset_name['test']:
            x_test_dnn_pred[name_test] = dnn_pred[name_test]

        rr_pred = np.load(os.path.join(args.inter_dir, 'rr_prediction.npz'), allow_pickle=True)
        x_train_rr_pred = rr_pred[src_dataset]
        x_test_rr_pred = {}
        for name_test in args.dataset_name['test']:
            x_test_rr_pred[name_test] = rr_pred[name_test]

        x_train = np.concatenate((x_train_dnn_pred, x_train_rr_pred), axis=1)
        x_test = {}
        for name_test in args.dataset_name['test']:
            x_test[name_test] = np.concatenate((x_test_dnn_pred[name_test], x_test_rr_pred[name_test]), axis=1)

        # prediction from RR models trained on large source datasets
        x_test_large_pred = {}
        for large_dataset in args.dataset_name['large']:
            n_phe_large = get_phe_num(args.in_dir, large_dataset)
            x_train_large_pred = np.zeros((x_train_dnn_pred.shape[0], n_phe_large))
            for name_test in args.dataset_name['test']:
                x_test_large_pred[name_test] = np.zeros((x_test_dnn_pred[name_test].shape[0], n_phe_large))
                for i in range(n_phe_large):
                    large_pred = np.load(
                        os.path.join(args.inter_dir, large_dataset + '_rr_base_' + str(i) + '.npz'),
                        allow_pickle=True)
                    x_train_large_pred[:, i] = large_pred[src_dataset]
                    x_test_large_pred[name_test][:, i] = large_pred[name_test]

                x_test[name_test] = np.concatenate((x_test[name_test], x_test_large_pred[name_test]), axis=1)
            x_train = np.concatenate((x_train, x_train_large_pred), axis=1)
        npz = np.load(os.path.join(args.in_dir, src_dataset, src_dataset + '_dnn_input.npz'), allow_pickle=True)
        y_train_raw = npz['y_raw']
        y_train = np.array(y_train_raw, dtype=np.float64)

    y_train = y_train[:, phe_idx]

    nan_idx = np.isnan(y_train)
    x_train = x_train[~nan_idx, :]
    y_train = y_train[~nan_idx]
    # print(x_train.shape, x_test['HCP'].shape, x_test['HCPA'].shape)

    # z normalize based on y_train, t_sigma is the standard deviation of non-NaN value
    y_train, _, _, _ = misc_z_norm(y_train)

    _, test_pred = RR(x_train, y_train, None, x_test, args)

    # save results for future use
    if args.multilayer:
        file_str = src_dataset + '_rr_multilayer_' + str(phe_idx) + '.npz'
    else:
        file_str = src_dataset + '_rr_base_' + str(phe_idx) + '.npz'
    result_str = os.path.join(args.inter_dir, file_str)
    results = test_pred
    np.savez(result_str, **results)


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
    parser.add_argument('--model_dir', type=str, default=config.MODEL_DIR)
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--phe_idx', type=int, default=0)
    parser.add_argument('--src_dataset', type=str, default='GSP')
    parser.add_argument('--dataset_name', type=dict, default=config.DATASET_NAME)

    multilayer_parser = parser.add_mutually_exclusive_group(required=False)
    multilayer_parser.add_argument('--2layer', dest='multilayer', action='store_true')
    multilayer_parser.add_argument('--1layer', dest='multilayer', action='store_false')
    parser.set_defaults(multilayer=False)

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
