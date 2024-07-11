#!/usr/bin/env python3from sklearn.kernel_ridge import KernelRidge
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import pickle
import random
import argparse
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV
from config import config
from CBIG_misc import misc_z_norm, demean_normalize


def LRR(x_train, y_train, args):
    '''training LRR on extra-large dataset (e.g., UK Biobank) and save models for later use

    Args:
        x_train: input training data (functional connectivity vectors of training participants)
        y_train: target training label (phenotype of training participants)
        args: arguments from command line

    Returns:
        None
    '''
    param_grid = {
        'alpha': [
            0.00001, 0.0001, 0.001, 0.004, 0.007, 0.01, 0.04, 0.07, 0.1, 0.4,
            0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 15, 20
        ]
    }
    rr = GridSearchCV(Ridge(), cv=config.N_CV, param_grid=param_grid)
    rr.fit(demean_normalize(x_train), y_train)

    pickle.dump(rr, open(os.path.join(args.model_dir, args.src_dataset + '_rr_model_'
                                      + str(args.phe_idx) + '.sav'), 'wb'))

    return


def train(args):
    '''training KRR base models on UK Biobank dataset

    Args:
        args: args from command line

    Returns:
        None
    '''
    # set all the seed
    random.seed(args.seed)
    np.random.seed(args.seed)
    phe_idx = args.phe_idx

    if args.exp_dataset:
        args.in_dir = config.IN_DIR_UT if args.unit_test else config.IN_DIR_EXP
        args.out_dir = config.OUT_DIR_UT if args.unit_test else config.OUT_DIR_EXP
        args.inter_dir = config.INTER_DIR_UT if args.unit_test else config.INTER_DIR_EXP
        args.model_dir = config.MODEL_DIR_UT if args.unit_test else config.MODEL_DIR_EXP
        args.src_dataset = config.DATASET_NAME_EXP['extra-large']

    os.makedirs(args.model_dir, exist_ok=True)

    npz = np.load(os.path.join(args.in_dir, args.src_dataset, args.src_dataset + '_dnn_input.npz'), allow_pickle=True)
    x_train = npz['x_raw']
    y_train = npz['y_raw']
    y_train = np.array(y_train, dtype=np.float64)

    # train a KRR for each phenotype
    y_train = y_train[:, phe_idx]

    nan_idx = np.isnan(y_train)
    x_train = x_train[~nan_idx, :]
    y_train = y_train[~nan_idx]

    # z normalize based on y_train, t_sigma is the standard deviation of non-NaN value
    y_train, _, _, _ = misc_z_norm(y_train)

    LRR(x_train, y_train, args)


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--in_dir', type=str, default=config.IN_DIR)
    parser.add_argument('--out_dir', '-o', type=str, default=config.OUT_DIR)
    parser.add_argument('--inter_dir', type=str, default=config.INTER_DIR)
    parser.add_argument('--model_dir', type=str, default=config.MODEL_DIR)
    parser.add_argument('--src_dataset', type=str, default='UKBB')
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--phe_idx', type=int, default=0)

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
    train(get_args())
