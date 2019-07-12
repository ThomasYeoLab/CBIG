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
from pathlib import PurePath

import torch
import torch.utils.data
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from torch.optim.lr_scheduler import MultiStepLR

from config import config
from CBIG_model_pytorch import bcn, CBIG_dataset, bcn_tmp
from CBIG_mics import mics_infer_metric, mics_log


def train(args):
    '''main function for FNN network

    Args:
        args: args from command line

    Returns:
        None
    '''

    t_overall = time.time()
    print('\nCBIG BrainNetCNN for UKBB and sex with argument: ' + str(args))

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.enabled = False

    # set gpu number
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # load data
    npz = PurePath(args.path_data, 'data_brainnetcnn_sex.npz').as_posix()
    npz = np.load(npz)
    train_x = npz['train_x']
    train_y = npz['train_y']
    valid_x = npz['valid_x']
    valid_y = npz['valid_y']
    test_x = npz['test_x']
    test_y = npz['test_y']

    # load dataset for PyTorch
    dset_train = CBIG_dataset(train_x, train_y, for_sex=True)
    trainloader = DataLoader(
        dset_train, batch_size=args.batch_size, shuffle=True, num_workers=8)
    dset_valid = CBIG_dataset(valid_x, valid_y, for_sex=True)
    validLoader = DataLoader(
        dset_valid, batch_size=args.batch_size, shuffle=True, num_workers=8)
    dset_test = CBIG_dataset(test_x, test_y, for_sex=True)
    testLoader = DataLoader(
        dset_test, batch_size=args.batch_size, shuffle=False, num_workers=8)

    runs = args.runs  # numbers of ensemble runs
    epochs = args.epochs  # numbers of epochs per run

    # initialization of result record
    tra_los_record = np.zeros((runs, epochs))
    val_los_record = np.zeros((runs, epochs))
    tes_los_record = np.zeros((runs, epochs))
    val_auc_record = np.zeros((runs, epochs))
    tes_auc_record = np.zeros((runs, epochs))
    tra_auc_record = np.zeros((runs, epochs))
    tes_res_record = np.zeros((runs, epochs, test_x.shape[0], 2))
    final_original = None

    # Code running - with multiple ensemble runs
    for run in range(runs):

        # initialization of network
        if args.alpha is None:
            net = bcn(
                args.e2e,
                args.e2n,
                args.n2g,
                args.dropout,
                args.dim,
                for_sex=True)
        else:
            net = bcn_tmp(
                args.e2e,
                args.e2n,
                args.n2g,
                args.dropout,
                args.dim,
                args.alpha,
                for_sex=True)
        net.to(device)

        # other components of network
        criterion = nn.CrossEntropyLoss()
        optimizer = optim.SGD(
            net.parameters(),
            lr=args.lr,
            momentum=args.momentum,
            weight_decay=args.weight_decay)
        scheduler = MultiStepLR(
            optimizer,
            milestones=[args.scheduler_decrease, args.scheduler_decrease * 2],
            gamma=0.1)

        # start epoch training
        for epoch in range(epochs):
            scheduler.step()

            # training
            train_loss = 0.0
            net.train(True)
            for (x, y) in trainloader:
                x, y = x.to(device), y.to(device)
                optimizer.zero_grad()
                outputs = net(x)
                loss = criterion(outputs, y)
                loss.backward()
                optimizer.step()
                train_loss += loss.item()
            tra_los_record[run, epoch] = train_loss / len(trainloader)

            net.train(False)
            # Training
            auc, _ = mics_infer_metric(trainloader, net, criterion, device)
            tra_auc_record[run, epoch] = auc

            # validation
            auc, loss = mics_infer_metric(validLoader, net, criterion, device)
            val_auc_record[run, epoch] = auc
            val_los_record[run, epoch] = loss

            # test
            auc, loss, real, pred = mics_infer_metric(
                testLoader, net, criterion, device, need_value=True)
            if final_original is not None:
                assert np.array_equal(final_original, real)
            else:
                final_original = real
            tes_res_record[run, epoch, :, :] = pred
            tes_auc_record[run, epoch] = auc
            tes_los_record[run, epoch] = loss

    log_args = {
        'tra_los_record': tra_los_record,
        'val_los_record': val_los_record,
        'tes_los_record': tes_los_record,
        'tra_auc_record': tra_auc_record,
        'val_auc_record': val_auc_record,
        'tes_auc_record': tes_auc_record,
        'tes_res_record': tes_res_record,
        'final_original': final_original
    }

    mics_log('brainnetcnn', args.out_path, index=args.index, **log_args)
    print("time spent: {:.4f}".format(time.time() - t_overall))

    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()
    # general parameters
    parser.add_argument('--path_data', type=str, default=config.UKBB_INTER_DIR)
    parser.add_argument('--out_path', '-o', type=str, default=config.OUT_PATH)
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument(
        '--batch_size', type=int, default=config.UKBB_BATCH_SIZE)
    parser.add_argument('--epochs', type=int, default=config.UKBB_EPOCHS)
    parser.add_argument('--runs', type=int, default=config.UKBB_RUNS)
    parser.add_argument('--dim', type=int, default=config.UKBB_N_DIM)
    parser.add_argument('--gpu', type=int, default=0)

    # hyperparameter
    parser.add_argument('--index', type=int, default=None)
    parser.add_argument('--lr', type=float, default=1e-2)
    parser.add_argument('--momentum', type=float, default=0.9)
    parser.add_argument('--weight_decay', type=float, default=0.02)
    parser.add_argument('--scheduler_decrease', type=int, default=75)
    parser.add_argument('--dropout', type=float, default=0.5)
    parser.add_argument('--alpha', type=float, default=None)
    parser.add_argument('--e2e', type=int, default=8)
    parser.add_argument('--e2n', type=int, default=16)
    parser.add_argument('--n2g', type=int, default=64)

    return parser.parse_args()


if __name__ == '__main__':
    train(get_args())
