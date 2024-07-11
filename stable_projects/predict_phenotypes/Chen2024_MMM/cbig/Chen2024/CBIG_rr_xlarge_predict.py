#!/usr/bin/env python3
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
from config import config
from CBIG_misc import demean_normalize, get_phe_num


def predict(args):
    '''use LRR models (trained on UK Biobank) to predict source phenotypes for other datasets

    Args:
        args: args from command line

    Returns:
        None
    '''

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)

    dataset_name = config.DATASET_NAME

    if args.exp_dataset:
        args.in_dir = config.IN_DIR_UT if args.unit_test else config.IN_DIR_EXP
        args.out_dir = config.OUT_DIR_UT if args.unit_test else config.OUT_DIR_EXP
        args.inter_dir = config.INTER_DIR_UT if args.unit_test else config.INTER_DIR_EXP
        args.model_dir = config.MODEL_DIR_UT if args.unit_test else config.MODEL_DIR_EXP
        args.src_dataset = config.DATASET_NAME_EXP['extra-large']
        dataset_name = config.DATASET_NAME_EXP

    x_dict = {}
    pred_dict = {}
    n_phe_src = get_phe_num(args.in_dir, args.src_dataset)  # number of source phenotypes in UK Biobank
    for dataset in dataset_name['large'] + dataset_name['medium'] + dataset_name['test']:
        x_dict[dataset] = np.load(os.path.join(args.in_dir, dataset, dataset + '_dnn_input.npz'))['x_raw']
        x_dict[dataset][np.isnan(x_dict[dataset])] = 0
        pred_dict[dataset] = np.zeros((x_dict[dataset].shape[0], n_phe_src))

    for phe_idx in range(n_phe_src):
        kr = pickle.load(open(os.path.join(args.model_dir,
            args.src_dataset + '_rr_model_' + str(phe_idx) + '.sav'), 'rb'))
        for dataset in dataset_name['large'] + dataset_name['medium'] + dataset_name['test']:
            pred_dict[dataset][:, phe_idx] = kr.predict(demean_normalize(x_dict[dataset]))

    file_str = 'rr_prediction.npz'
    result_str = os.path.join(args.inter_dir, file_str)
    np.savez(result_str, **pred_dict)


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--in_dir', type=str, default=config.IN_DIR)
    parser.add_argument('--out_dir', '-o', type=str, default=config.OUT_DIR)
    parser.add_argument('--inter_dir', type=str, default=config.INTER_DIR)
    parser.add_argument('--model_dir', type=str, default=config.MODEL_DIR)
    parser.add_argument('--src_dataset', type=str, default='UKBB')
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)

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
    predict(get_args())
