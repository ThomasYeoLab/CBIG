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
import torch
import numpy as np
import torch.utils.data
from torch.utils.data import DataLoader
from config import config
from CBIG_misc import demean_normalize


def predict(args):
    '''main function for DNN network

    Args:
        args: args from command line

    Returns:
        None
    '''

    t_overall = time.time()

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    if args.exp_dataset:
        args.in_dir = config.IN_DIR_UT if args.unit_test else config.IN_DIR_EXP
        args.out_dir = config.OUT_DIR_UT if args.unit_test else config.OUT_DIR_EXP
        args.inter_dir = config.INTER_DIR_UT if args.unit_test else config.INTER_DIR_EXP
        args.model_dir = config.MODEL_DIR_UT if args.unit_test else config.MODEL_DIR_EXP
        args.src_dataset = config.DATASET_NAME_EXP['extra-large']

    # set gpu number
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # load data
    npz = np.load(os.path.join(args.inter_dir, 'dnn_base.npz'))
    n_phe = npz['tra_cor_record'].shape[1]
    tuned_by = args.tuned_by
    val_record = npz['val_' + tuned_by + '_record']
    temp = np.mean(val_record[:, :], axis=1)
    temp = np.convolve(temp, np.ones(3, dtype=int), 'valid') / 3
    index = np.nanargmax(temp)
    index = index + 1
    print('\nBest validation at index: ', index)
    print(np.mean(val_record[index, :]))

    # load model
    dir_model = os.path.join(args.model_dir, 'trained_model_' + args.src_dataset)
    dir_model_temp = os.path.join(dir_model, 'dnn_model_save_base')
    model_path = os.path.join(dir_model_temp, 'CBIG_dnn_epoch_' + str(index) + '.pkl_torch')
    net = torch.load(model_path)
    net.to(device)
    net.train(False)

    dataset_name = config.DATASET_NAME_EXP if args.exp_dataset else config.DATASET_NAME
    x_dict = {}
    loader_dict = {}
    res_record_dict = {}

    for dataset in dataset_name['large'] + dataset_name['medium'] + dataset_name['test']:
        x_dict[dataset] = np.load(os.path.join(args.in_dir, dataset, dataset + '_dnn_input.npz'))['x_raw']
        x_dict[dataset][np.isnan(x_dict[dataset])] = 0
        # subject-wise normalization for input functional connectivity
        x_dict[dataset] = demean_normalize(x_dict[dataset]).astype(np.float32)
        # load dataset for PyTorch
        loader_dict[dataset] = DataLoader(x_dict[dataset], batch_size=args.batch_size, shuffle=False, num_workers=0)
        # initialization of result record
        res_record_dict[dataset] = np.zeros((x_dict[dataset].shape[0], n_phe))
        record_pred = np.zeros((0, n_phe))  # prediction value
        for x in loader_dict[dataset]:
            x = x.to(device)
            outputs = net(x)
            record_pred = np.concatenate((record_pred, outputs.data.cpu().numpy()), axis=0)
        res_record_dict[dataset] = np.squeeze(record_pred)

    model_str = 'dnn_prediction'
    np.savez(os.path.join(args.inter_dir, model_str + '.npz'), **res_record_dict)

    print("time spent: {:.4f}".format(time.time() - t_overall))


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
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--batch_size', type=int, default=config.BATCH_SIZE)
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--tuned_by', type=str, default='cod')
    parser.add_argument('--src_dataset', type=str, default='UKBB')

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
