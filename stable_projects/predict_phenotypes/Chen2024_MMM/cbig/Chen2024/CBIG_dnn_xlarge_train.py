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
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split
from config import config
from CBIG_model_pytorch import dnn, msenanloss, multi_task_dataset
from CBIG_misc import misc_z_norm, misc_infer_metric, misc_log, print_result, demean_normalize


def train(args):
    '''train DNN network on extra-large dataset (e.g., UK Biobank)

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

    print('\nCBIG DNN for extra-large dataset with argument: ' + str(args))


    # set gpu number
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # load data
    npz = np.load(os.path.join(args.in_dir, args.src_dataset, args.src_dataset + '_dnn_input.npz'))
    x_train_raw = npz['x_raw']
    x_train_raw[np.isnan(x_train_raw)] = 0
    y_train_raw = npz['y_raw']

    # subject-wise normalization for input functional connectivity
    x_train_raw = demean_normalize(x_train_raw)

    # split train and validation
    split_tra, split_val = train_test_split(np.arange(x_train_raw.shape[0]), test_size=0.2, random_state=seed)

    x_train = x_train_raw[split_tra, :]
    x_valid = x_train_raw[split_val, :]
    y_train = y_train_raw[split_tra, :]
    y_valid = y_train_raw[split_val, :]

    # z norm based on y_train
    y_train, y_valid, _, t_sigma = misc_z_norm(y_train, y_valid)

    # since we perform z normalization before, and we not need MAE now
    # so we use ones for t_sigma that we might need it in the future
    n_phe = y_train.shape[1]

    # load dataset for PyTorch
    dset_train = multi_task_dataset(x_train, y_train)
    trainloader = DataLoader(dset_train, batch_size=args.batch_size, shuffle=True, num_workers=0)
    dset_valid = multi_task_dataset(x_valid, y_valid)
    validLoader = DataLoader(dset_valid, batch_size=args.batch_size, shuffle=True, num_workers=0)

    epochs = args.epochs  # numbers of epochs

    # initialization of result record
    tra_los_record = np.zeros((epochs))
    val_los_record = np.zeros((epochs))
    tra_cor_record = np.zeros((epochs, n_phe))
    val_cor_record = np.zeros((epochs, n_phe))
    tra_cod_record = np.zeros((epochs, n_phe))
    val_cod_record = np.zeros((epochs, n_phe))
    tra_mae_record = np.zeros((epochs, n_phe))
    val_mae_record = np.zeros((epochs, n_phe))
    val_res_record = np.zeros((epochs, x_valid.shape[0], n_phe))

    os.makedirs(args.inter_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    dir_model = os.path.join(args.model_dir, 'trained_model_' + args.src_dataset)
    os.makedirs(dir_model, exist_ok=True)
    dir_model_temp = os.path.join(dir_model, 'dnn_model_save_base')
    os.makedirs(dir_model_temp, exist_ok=True)

    # initialization of network
    net = dnn(
        x_train.shape[1],
        args.n_hidden_layer,
        args.n_l1,
        args.n_l2,
        args.n_l3,
        args.n_l4,
        args.dropout,
        output_size=n_phe)
    print(net)
    net.to(device)

    # other components of network
    criterion = msenanloss
    optimizer = optim.SGD(
        net.parameters(),
        lr=args.lr,
        momentum=args.momentum,
        weight_decay=args.weight_decay)

    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=0.5,
        patience=args.patience,
        verbose=True)

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
        tra_los_record[epoch] = train_loss / len(trainloader)

        net.train(False)
        # Training
        corr, cod, mae, _ = misc_infer_metric(
            trainloader,
            net,
            criterion,
            device,
            t_sigma,
            output_size=n_phe)
        tra_cor_record[epoch, :] = corr
        tra_cod_record[epoch, :] = cod
        tra_mae_record[epoch, :] = mae

        # validation
        corr, cod, mae, loss, _, pred = misc_infer_metric(
            validLoader,
            net,
            criterion,
            device,
            t_sigma,
            output_size=n_phe,
            need_value=True)
        val_cor_record[epoch, :] = corr
        val_cod_record[epoch, :] = cod
        val_mae_record[epoch, :] = mae
        val_los_record[epoch] = loss
        val_res_record[epoch, :, :] = np.squeeze(pred)
        print_result(tra_cor_record[epoch, :], tra_cod_record[epoch, :],
                     tra_mae_record[epoch, :], loss, corr, cod, mae, epoch)

        # save model
        model_path = os.path.join(dir_model_temp, 'CBIG_dnn_epoch_' + str(epoch) + '.pkl_torch')
        torch.save(net, model_path)

        scheduler.step(loss)

    log_args = {
        'tra_los_record': tra_los_record,
        'val_los_record': val_los_record,
        'tra_cor_record': tra_cor_record,
        'val_cor_record': val_cor_record,
        'tra_cod_record': tra_cod_record,
        'val_cod_record': val_cod_record,
        'tra_mae_record': tra_mae_record,
        'val_mae_record': val_mae_record,
        't_sigma': t_sigma
    }

    model_str = 'dnn'
    misc_log(
        model_str,
        args.inter_dir,
        metric=args.metric,
        index=args.index,
        **log_args)
    print("time spent: {:.4f}".format(time.time() - t_overall))


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--in_dir', type=str, default=config.IN_DIR)
    parser.add_argument('--inter_dir', type=str, default=config.INTER_DIR)
    parser.add_argument('--out_dir', type=str, default=config.OUT_DIR)
    parser.add_argument('--model_dir', type=str, default=config.MODEL_DIR)
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--batch_size', type=int, default=config.BATCH_SIZE)
    parser.add_argument('--epochs', type=int, default=config.EPOCHS)
    parser.add_argument('--metric', type=str, default='cod')
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--index', type=int, default=None)
    parser.add_argument('--src_dataset', type=str, default='UKBB')

    # hyperparameter
    parser.add_argument('--lr', type=float, default=1e-2)
    parser.add_argument('--momentum', type=float, default=0.9)
    parser.add_argument('--weight_decay', type=float, default=1e-5)
    parser.add_argument('--patience', type=int, default=25)
    parser.add_argument('--dropout', type=float, default=0.2)
    parser.add_argument('--n_hidden_layer', type=int, default=3)
    parser.add_argument('--n_l1', type=int, default=64)
    parser.add_argument('--n_l2', type=int, default=32)
    parser.add_argument('--n_l3', type=int, default=32)
    parser.add_argument('--n_l4', type=int, default=32)

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
