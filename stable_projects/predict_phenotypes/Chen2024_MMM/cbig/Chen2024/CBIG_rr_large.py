#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import random
import pickle
import argparse
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import GridSearchCV
from config import config
from CBIG_misc import misc_z_norm, demean_normalize


def RR(x_train, y_train, x_medium, x_test, args):
    '''training RR on medium (or large) source dataset and predict on smaller datasets

    Args:
        x_train: input training data (functional connectivity vectors of training participants)
        y_train: target training label (phenotype of training participants)
        x_medium: functional connectivity vectors of participants on medium datasets (e.g., HBN, GSP, eNKI)
            x_medium = {name_medium: x_medium, ...}
        x_test: functional connectivity vectors of participants on meta-test datasets (e.g., HCP, HCP-Aging)
            x_test = {name_test: x_test, ...}
        args: arguments from command line

    Returns:
        medium_pred: prediction of phenotypes (from large source datasets) on medium source datasets
            medium_pred = {name_medium: medium_pred}
        test_pred: prediction of phenotypes (from large/medium source datasets) on meta-test datasets
            test_pred = {name_test: test_pred, ...}
    '''
    param_grid = {
        'alpha': [
            0.00001, 0.0001, 0.001, 0.004, 0.007, 0.01, 0.04, 0.07, 0.1, 0.4,
            0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 15, 20
        ]
    }
    if args.multilayer:
        model = KernelRidge()
        log_str = 'multilayer'
    else:
        model = Ridge()
        log_str = 'base'
    rr = GridSearchCV(model, cv=config.N_CV, param_grid=param_grid)
    rr.fit(demean_normalize(x_train), y_train)

    os.makedirs(args.model_dir, exist_ok=True)
    pickle.dump(rr, open(os.path.join(args.model_dir, args.src_dataset + '_rr_model_' + log_str + '_'
                                      + str(args.phe_idx) + '.sav'), 'wb'))

    test_pred = {}
    for name_test in args.dataset_name['test']:
        test_pred[name_test] = rr.predict(demean_normalize(x_test[name_test]))

    medium_pred = {}
    if x_medium != None:
        for name_medium in x_medium:
            medium_pred[name_medium] = rr.predict(demean_normalize(x_medium[name_medium]))

    return medium_pred, test_pred


def main(args):
    '''train Ridge Regression base models on large source dataset (e.g., ABCD) and predict on test datasets

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
        args.src_dataset = config.DATASET_NAME_EXP['large'][0]
        args.dataset_name = config.DATASET_NAME_EXP

    phe_idx = args.phe_idx
    src_dataset = args.src_dataset

    # load input and label data
    if args.multilayer == False:
        # apply direct meta-matching
        npz = np.load(os.path.join(args.in_dir, src_dataset, src_dataset + '_dnn_input.npz'), allow_pickle=True)
        x_train = npz['x_raw']
        x_train[np.isnan(x_train)] = 0
        y_train= npz['y_raw'][:, phe_idx]
        y_train = np.array(y_train, dtype=np.float64)
        x_test = {}  # read FC data in all test datasets
        for tar_dataset in args.dataset_name['test']:
            x_test[tar_dataset] = np.load(os.path.join(args.in_dir, tar_dataset, tar_dataset +
                                                         '_dnn_input.npz'), allow_pickle=True)['x_raw']
            x_test[tar_dataset][np.isnan(x_test[tar_dataset])] = 0

        x_medium = {} # read FC data in all medium source datasets
        for name_medium in args.dataset_name['medium']:
            x_medium[name_medium] = np.load(os.path.join(args.in_dir, name_medium, name_medium + '_dnn_input.npz'),
                                            allow_pickle=True)['x_raw']
            x_medium[name_medium][np.isnan(x_medium[name_medium])] = 0
    else:
        # perform multi-layer meta-matching (2-layer)
        dnn_pred = np.load(os.path.join(args.inter_dir, 'dnn_prediction.npz'), allow_pickle=True)
        x_train_dnn_pred = dnn_pred[src_dataset]
        x_test_dnn_pred = {}
        for name_test in args.dataset_name['test']:
            x_test_dnn_pred[name_test] = dnn_pred[name_test]
        x_medium_dnn_pred = {}
        for name_medium in args.dataset_name['medium']:
            x_medium_dnn_pred[name_medium]= dnn_pred[name_medium]

        rr_pred = np.load(os.path.join(args.inter_dir, 'rr_prediction.npz'), allow_pickle=True)
        x_train_rr_pred = rr_pred[src_dataset]
        # combine RR prediction and DNN prediction
        x_train = np.concatenate((x_train_dnn_pred, x_train_rr_pred), axis=1)
        x_test_rr_pred = {}
        x_test = {}  # combine RR prediction and DNN prediction
        for name_test in args.dataset_name['test']:
            x_test_rr_pred[name_test] = rr_pred[name_test]
            x_test[name_test] = np.concatenate((x_test_dnn_pred[name_test], x_test_rr_pred[name_test]), axis=1)
        x_medium_rr_pred = {}
        x_medium = {} # combine RR prediction and DNN prediction
        for name_medium in args.dataset_name['medium']:
            x_medium_rr_pred[name_medium] = rr_pred[name_medium]
            x_medium[name_medium] = np.concatenate((x_medium_dnn_pred[name_medium],
                                                    x_medium_rr_pred[name_medium]), axis = 1)

        npz = np.load(os.path.join(args.in_dir, src_dataset, src_dataset + '_dnn_input.npz'), allow_pickle=True)
        y_train = npz['y_raw'][:, phe_idx]
        y_train = np.array(y_train, dtype=np.float64)

    nan_idx = np.isnan(y_train)
    x_train = x_train[~nan_idx, :]
    y_train = y_train[~nan_idx]

    # z normalize based on y_train
    y_train, _, _, _ = misc_z_norm(y_train)

    medium_pred, test_pred = RR(x_train, y_train, x_medium, x_test, args)

    # save results for future use
    if args.multilayer:
        file_str = src_dataset + '_rr_multilayer_' + str(phe_idx) + '.npz'
    else:
        file_str = src_dataset + '_rr_base_' + str(phe_idx) + '.npz'

    result_str = os.path.join(args.inter_dir, file_str)
    results = {**medium_pred, **test_pred}
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
    parser.add_argument('--src_dataset', type=str, default='ABCD')
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
