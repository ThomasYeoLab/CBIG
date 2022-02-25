#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import random
import argparse
import numpy as np
from sklearn.model_selection import train_test_split

import torch
import torch.utils.data
from torch.utils.data import DataLoader

from config import config
from CBIG_model_pytorch import ukbb_multi_task_dataset


def pred_y_train(x_train, y_train_dummy, index, device, args):
    '''pred y with training meta-set data

    Args:
        x_train (ndarray): training set of training meta-set FC data
        y_train_dummy (ndarray): dummy y data
        index (int): optimal dnn epoch index for base training
        device (torch.device): the gpu that torch runs on
        args: args from command line

    Returns:
        ndarray: predicted y for haufe transform
    '''
    dset_train = ukbb_multi_task_dataset(x_train, y_train_dummy, True)
    trainloader = DataLoader(
        dset_train, batch_size=128, shuffle=False, num_workers=8)
    weight_path = os.path.join(
        args.out_dir, 'trained_model_ukbb',
        'dnn_model_save_base_cross_dataset',
        'CBIG_ukbb_dnn_run_0_epoch_' + str(index) + '.pkl_torch')
    net = torch.load(weight_path)

    net.train(False)
    record_pred = np.zeros((0, y_train_dummy.shape[1]))
    for (x, _) in trainloader:
        x = x.to(device)
        outputs = net(x)
        record_pred = np.concatenate((record_pred, outputs.data.cpu().numpy()),
                                     axis=0)
    return np.squeeze(record_pred)


def get_haufe_ukbb_training(args):
    '''Get data for haufe transform on UK Biobank training meta-set

    Args:
        args: args from command line

    Returns:
        None
    '''

    print('\nArgument: ' + str(args))

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    tuned_by = 'cod'
    os.makedirs(os.path.join(args.out_dir, 'tmp'), exist_ok=True)

    # load base prediction result
    npz = os.path.join(args.out_dir, 'dnn_across_dataset_base.npz')
    npz = np.load(npz)
    val_record = npz['val_' + tuned_by + '_record']
    temp = np.mean(val_record[0, :, :], axis=1)
    temp = np.convolve(temp, np.ones(3, dtype=int), 'valid') / 3
    index = np.nanargmax(temp)
    index = index + 1
    print('\nBest validation at index: ', index)

    # load original data
    npz = os.path.join(args.in_dir, 'ukbb_dnn_input_cross_dataset.npz')
    npz = np.load(npz, allow_pickle=True)
    tra_phe = npz['tra_phe']
    x_train_raw = npz['x_train_raw']
    y_train_raw = npz['y_train_raw']
    split_tra, _ = train_test_split(
        np.arange(x_train_raw.shape[0]), test_size=0.2, random_state=seed)
    x_train = x_train_raw[split_tra, :]
    y_train = y_train_raw[split_tra, :]
    y_train_dummy = np.zeros(y_train.shape)

    y_pred_train = pred_y_train(x_train, y_train_dummy, index, device, args)
    npz_train_pred = os.path.join(args.out_dir, 'haufe_y_pred_train.npz')
    np.savez(
        npz_train_pred,
        y_pred_train=y_pred_train,
        x_train=x_train,
        tra_phe=tra_phe)

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
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--split', type=str, default='test')

    return parser.parse_args()


if __name__ == '__main__':
    get_haufe_ukbb_training(get_args())
