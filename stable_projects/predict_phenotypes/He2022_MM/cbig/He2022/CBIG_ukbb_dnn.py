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
import scipy.io as sio
from sklearn.model_selection import train_test_split

import torch
import torch.utils.data
import torch.optim as optim
from torch.utils.data import DataLoader
from torch.optim.lr_scheduler import MultiStepLR

from config import config
from CBIG_model_pytorch import dnn_4l, dnn_3l, dnn_2l, dnn_5l
from CBIG_model_pytorch import msenanloss, ukbb_multi_task_dataset
from CBIG_mics import mics_z_norm, mics_infer_metric, mics_log, print_result


def train(args):
    '''main function for DNN network

    Args:
        args: args from command line

    Returns:
        None
    '''

    t_overall = time.time()
    print('\nCBIG DNN for UK Biobank with argument: ' + str(args))

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    # set gpu number
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # load data
    if args.across_dataset:
        npz = os.path.join(args.in_dir, 'ukbb_dnn_input_cross_dataset.npz')
    else:
        npz = os.path.join(args.in_dir, 'ukbb_dnn_input_test.npz')
    npz = np.load(npz)
    x_train_raw = npz['x_train_raw']
    y_train_raw = npz['y_train_raw']
    x_test = npz['x_test']
    y_test_different_set_phe = npz['y_test_different_set_phe']

    # split train and validation
    if args.across_dataset:
        split_tra, split_val = train_test_split(
            np.arange(x_train_raw.shape[0]), test_size=0.2, random_state=seed)
    else:
        split_file = os.path.join(args.in_dir,
                                  'split_rng' + str(seed) + '.mat')
        split = np.squeeze(sio.loadmat(split_file)['sub_fold'][0][0][0])
        split_tra = split == 0
        split_val = split == 2
        # split_tes = split == 1
    x_train = x_train_raw[split_tra, :]
    x_valid = x_train_raw[split_val, :]
    y_train = y_train_raw[split_tra, :]
    y_valid = y_train_raw[split_val, :]

    # z norm based on y_train
    y_train, y_valid, _, t_sigma = mics_z_norm(y_train, y_valid)

    # since we perform z normalization before, and we not need MAE now
    # so we use ones for t_sigma that we might need it in the future
    n_phe = y_train.shape[1]

    # load dataset for PyTorch
    dset_train = ukbb_multi_task_dataset(x_train, y_train)
    trainloader = DataLoader(
        dset_train, batch_size=args.batch_size, shuffle=True, num_workers=8)
    dset_valid = ukbb_multi_task_dataset(x_valid, y_valid)
    validLoader = DataLoader(
        dset_valid, batch_size=args.batch_size, shuffle=True, num_workers=8)
    dset_test = ukbb_multi_task_dataset(x_test, y_test_different_set_phe)
    testLoader = DataLoader(
        dset_test, batch_size=args.batch_size, shuffle=False, num_workers=8)

    runs = args.runs  # numbers of ensemble runs
    epochs = args.epochs  # numbers of epochs per run

    # initialization of result record
    tra_los_record = np.zeros((runs, epochs))
    val_los_record = np.zeros((runs, epochs))
    tra_cor_record = np.zeros((runs, epochs, n_phe))
    val_cor_record = np.zeros((runs, epochs, n_phe))
    tra_cod_record = np.zeros((runs, epochs, n_phe))
    val_cod_record = np.zeros((runs, epochs, n_phe))
    tra_mae_record = np.zeros((runs, epochs, n_phe))
    val_mae_record = np.zeros((runs, epochs, n_phe))
    tes_res_record = np.zeros((runs, epochs, x_test.shape[0], n_phe))
    val_res_record = np.zeros((runs, epochs, x_valid.shape[0], n_phe))
    final_original = None

    os.makedirs(args.out_dir, exist_ok=True)
    dir_model = os.path.join(args.out_dir, 'trained_model_ukbb')
    os.makedirs(dir_model, exist_ok=True)
    dir_model_temp = os.path.join(dir_model, 'dnn_model_save_base')
    if args.across_dataset:
        dir_model_temp += '_cross_dataset'
    os.makedirs(dir_model_temp, exist_ok=True)

    # Code running - with multiple ensemble runs
    for run in range(runs):

        # initialization of network
        if args.n_layer == 2:
            net = dnn_2l(
                x_train.shape[1], args.n_l1, args.dropout, output_size=n_phe)
        elif args.n_layer == 3:
            net = dnn_3l(
                x_train.shape[1],
                args.n_l1,
                args.n_l2,
                args.dropout,
                output_size=n_phe)
        elif args.n_layer == 4:
            net = dnn_4l(
                x_train.shape[1],
                args.n_l1,
                args.n_l2,
                args.n_l3,
                args.dropout,
                output_size=n_phe)
        elif args.n_layer == 5:
            net = dnn_5l(
                x_train.shape[1],
                args.n_l1,
                args.n_l2,
                args.n_l3,
                args.n_l4,
                args.dropout,
                output_size=n_phe)
        else:
            assert False, "Only support 2, 3, 4, 5 layers."
        net.to(device)

        # other components of network
        criterion = msenanloss
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

            # training
            train_loss = 0.0
            net.train(True)
            for (x, y) in trainloader:
                x, y = x.to(device), y.to(device)
                optimizer.zero_grad()
                outputs = net(x)
                mask = torch.isnan(y)
                y.masked_fill_(mask, 0)
                loss = criterion(outputs, y, mask=mask)
                loss.backward()
                optimizer.step()
                train_loss += loss.item()
            tra_los_record[run, epoch] = train_loss / len(trainloader)

            net.train(False)
            # Training
            corr, cod, mae, _ = mics_infer_metric(
                trainloader,
                net,
                criterion,
                device,
                t_sigma,
                output_size=n_phe)
            tra_cor_record[run, epoch, :] = corr
            tra_cod_record[run, epoch, :] = cod
            tra_mae_record[run, epoch, :] = mae

            # validation
            corr, cod, mae, loss, _, pred = mics_infer_metric(
                validLoader,
                net,
                criterion,
                device,
                t_sigma,
                output_size=n_phe,
                need_value=True)
            val_cor_record[run, epoch, :] = corr
            val_cod_record[run, epoch, :] = cod
            val_mae_record[run, epoch, :] = mae
            val_los_record[run, epoch] = loss
            val_res_record[run, epoch, :, :] = np.squeeze(pred)
            print_result(
                tra_cor_record[run, epoch, :], tra_cod_record[run, epoch, :],
                tra_mae_record[run, epoch, :], loss, corr, cod, mae, epoch)

            record_pred = np.zeros((0, n_phe))  # prediction value
            for (x, y) in testLoader:
                x, y = x.to(device), y.to(device)
                outputs = net(x)
                record_pred = np.concatenate(
                    (record_pred, outputs.data.cpu().numpy()), axis=0)
            tes_res_record[run, epoch, :, :] = np.squeeze(record_pred)

            file_trained = os.path.join(
                dir_model_temp, 'CBIG_ukbb_dnn_run_' + str(run) + '_epoch_' +
                str(epoch) + '.pkl_torch')
            torch.save(net, file_trained)

            scheduler.step()

    log_args = {
        'tra_los_record': tra_los_record,
        'val_los_record': val_los_record,
        'tra_cor_record': tra_cor_record,
        'val_cor_record': val_cor_record,
        'tra_cod_record': tra_cod_record,
        'val_cod_record': val_cod_record,
        'tra_mae_record': tra_mae_record,
        'val_mae_record': val_mae_record,
        'tes_res_record': tes_res_record,
        'final_original': final_original,
        'val_res_record': val_res_record,
        't_sigma': t_sigma
    }

    model_str = 'dnn'
    if args.across_dataset:
        model_str += '_across_dataset'
    mics_log(
        model_str,
        args.out_dir,
        metric=args.metric,
        index=args.index,
        **log_args)
    print("time spent: {:.4f}".format(time.time() - t_overall))

    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--in_dir', type=str, default=config.IN_DIR)
    parser.add_argument('--out_dir', '-o', type=str, default=config.OUT_DIR)
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--batch_size', type=int, default=config.BATCH_SIZE)
    parser.add_argument('--epochs', type=int, default=config.EPOCHS)
    parser.add_argument('--runs', type=int, default=config.RUNS)
    parser.add_argument('--metric', type=str, default='cod')
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--across_dataset', type=bool, default=False)

    # hyperparameter
    parser.add_argument('--index', type=int, default=None)
    parser.add_argument('--lr', type=float, default=1e-2)
    parser.add_argument('--momentum', type=float, default=0.9)
    parser.add_argument('--weight_decay', type=float, default=0.02)
    parser.add_argument('--scheduler_decrease', type=int, default=75)
    parser.add_argument('--dropout', type=float, default=0.5)
    parser.add_argument('--n_layer', type=int, default=3)
    parser.add_argument('--n_l1', type=int, default=64)
    parser.add_argument('--n_l2', type=int, default=32)
    parser.add_argument('--n_l3', type=int, default=32)
    parser.add_argument('--n_l4', type=int, default=32)

    return parser.parse_args()


if __name__ == '__main__':
    train(get_args())
