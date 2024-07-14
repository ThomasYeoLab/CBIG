#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import time
import random
import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

import torch
import torch.utils.data
import torch.optim as optim
from torch.utils.data import DataLoader
from torch.optim.lr_scheduler import MultiStepLR

from cbig.CBIG_model_pytorch import SFCN, msenanloss, vol_dataset
from cbig.CBIG_mics import mics_z_norm, mics_infer_metric, \
    mics_log, save_model, read_datapath
from cbig.config import config


def train(args):
    '''main function for training base DNN network (this is
    basic meta-matching (DNN) for training)

    Args:
        args: args from command line

    Returns:
        None
    '''

    t_overall = time.time()  # seconds
    print('\nCBIG DNN for UK Biobank with argument: ' + str(args))

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.enabled = True

    if not os.path.exists(
            os.path.join(args.out_dir, args.out_subdir + '/training_process')):
        os.makedirs(
            os.path.join(args.out_dir, args.out_subdir + '/training_process'))

    # set gpu
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # load data
    phe = pd.read_csv(os.path.join(args.phe_dir))
    icv = pd.read_csv(os.path.join(args.icv_dir))
    sub_list = pd.read_table(os.path.join(args.sub_dir), header=None)
    mask = phe['eid'].isin(sub_list[0].values.tolist())
    phe = phe.loc[mask, phe.columns[args.start_idx:args.end_idx]]
    mask_icv = icv['eid'].isin(sub_list[0].values.tolist())
    icv = icv.loc[mask_icv, ['inverse_determinant']]
    y_train_raw = phe.values.astype(float)
    icv_raw = icv.values.astype(float)
    x_train_raw = read_datapath(args.data_dir, list(map(str, sub_list[0])),
                                args.dataset)

    # split train and validation
    split_tra, split_val = train_test_split(np.arange(x_train_raw.shape[0]),
                                            test_size=0.2,
                                            random_state=seed)
    x_train = x_train_raw[split_tra]
    x_valid = x_train_raw[split_val]
    y_train = y_train_raw[split_tra, :]
    y_valid = y_train_raw[split_val, :]
    icv_train = icv_raw[split_tra, :]
    icv_valid = icv_raw[split_val, :]
    np.savez(os.path.join(
        args.out_dir,
        args.out_subdir + '/training_process/train_val_split_ukbb.npz'),
             split_tra=split_tra,
             split_val=split_val)

    # z norm
    y_train, y_valid, _, t_sigma = mics_z_norm(y_train, y_valid)
    icv_train, icv_valid, _, _ = mics_z_norm(icv_train, icv_valid)

    n_phe = y_train.shape[1]

    # load dataset
    batch_size = args.batch_size
    dset_train = vol_dataset(x_train, y_train, icv=icv_train)
    trainloader = DataLoader(dset_train,
                             batch_size=batch_size,
                             shuffle=True,
                             num_workers=batch_size)
    dset_valid = vol_dataset(x_valid, y_valid, icv=icv_valid)
    validLoader = DataLoader(dset_valid,
                             batch_size=batch_size,
                             shuffle=True,
                             num_workers=batch_size)

    runs = args.runs
    epochs = args.epochs  # numbers of epochs per run
    pre_index = 0  # record epoch index for best model

    # initialization of result record
    tra_los_record = np.zeros((runs, epochs))
    val_los_record = np.zeros((runs, epochs))
    val_cor_record = np.zeros((runs, epochs, n_phe))
    val_cod_record = np.zeros((runs, epochs, n_phe))
    val_mae_record = np.zeros((runs, epochs, n_phe))
    val_res_record = np.zeros((runs, epochs, x_valid.shape[0], n_phe))

    os.makedirs(args.out_dir, exist_ok=True)
    dir_model = os.path.join(args.out_dir, args.out_subdir)
    os.makedirs(dir_model, exist_ok=True)
    dir_model_temp = os.path.join(dir_model, 'dnn_model_save_base')
    if args.across_dataset:
        dir_model_temp += '_cross_dataset'
    os.makedirs(dir_model_temp, exist_ok=True)

    for run in range(runs):

        net = SFCN(channel_number=args.channel_number,
                   output_dim=args.output_dim,
                   dropout=args.dropout)
        net.to(device)

        # other components of network
        criterion = msenanloss

        # original optimizer
        optimizer = optim.SGD(net.parameters(),
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
            for (x, y, icv) in trainloader:
                x, y, icv = x.to(device), y.to(device), icv.to(device)
                optimizer.zero_grad()
                outputs = net(x, icv)
                mask = torch.isnan(y)
                y.masked_fill_(mask, 0)
                loss = criterion(outputs, y, mask=mask)

                del outputs, x, y, icv
                torch.cuda.empty_cache()

                loss.backward()
                optimizer.step()
                train_loss += loss.item()

                del loss
                torch.cuda.empty_cache()
            tra_los_record[run, epoch] = train_loss / len(trainloader)

            net.train(False)
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

            scheduler.step()

            if epoch >= 1:
                file_trained = os.path.join(
                    dir_model_temp, 'CBIG_ukbb_dnn_run_' + str(run) +
                    '_epoch_' + str(epoch) + '.pkl_torch')
                val_log_args = {
                    'val_los_record': val_los_record,
                    'val_cor_record': val_cor_record,
                    'val_cod_record': val_cod_record
                }
                pre_index = save_model(net,
                                       file_trained,
                                       pre_index,
                                       metric=args.metric,
                                       **val_log_args)

    log_args = {
        'tra_los_record': tra_los_record,
        'val_los_record': val_los_record,
        'val_cor_record': val_cor_record,
        'val_cod_record': val_cod_record,
        'val_mae_record': val_mae_record,
        'val_res_record': val_res_record,
        't_sigma': t_sigma
    }

    model_str = 'dnn'
    if args.across_dataset:
        model_str += '_across_dataset'
    mics_log(model_str,
             os.path.join(args.out_dir, args.out_subdir + '/training_process'),
             metric=args.metric,
             index=args.index,
             **log_args)
    print("time spent: {:.4f}".format(time.time() - t_overall))

    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by
          other function

    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--in_dir', type=str, default=None)
    parser.add_argument('--phe_dir', type=str, default=None)
    parser.add_argument('--icv_dir', type=str, default=None)
    parser.add_argument('--sub_dir', type=str, default=None)
    parser.add_argument('--data_dir', type=str, default=None)
    parser.add_argument('--out_dir', '-o', type=str, default=None)
    parser.add_argument('--out_subdir', '-osub', type=str, default=None)
    parser.add_argument('--inter_dir', type=str, default=None)
    parser.add_argument('--seed', type=int, default=config.SEED)
    parser.add_argument('--batch_size', type=int, default=None)
    parser.add_argument('--epochs', type=int, default=None)
    parser.add_argument('--runs', type=int, default=config.RUNS)
    parser.add_argument('--metric', type=str, default='cod')
    across_dataset_parser = parser.add_mutually_exclusive_group(required=False)
    across_dataset_parser.add_argument('--across-dataset',
                                       dest='across_dataset',
                                       action='store_true')
    across_dataset_parser.add_argument('--not-across-dataset',
                                       dest='across_dataset',
                                       action='store_false')
    parser.set_defaults(across_dataset=False)

    # hyperparameter
    parser.add_argument('--index', type=int, default=None)
    parser.add_argument('--lr', type=float, default=None)
    parser.add_argument('--momentum', type=float, default=config.MOMENTUM)
    parser.add_argument('--weight_decay',
                        type=float,
                        default=config.WEIGHT_DECAY)
    parser.add_argument('--scheduler_decrease',
                        type=int,
                        default=config.SCHEDULER_DECREASE)
    parser.add_argument("--channel_number",
                        "-cn",
                        type=int,
                        nargs='+',
                        default=None)
    parser.add_argument('--output_dim', type=int, default=None)
    parser.add_argument('--dropout', type=float, default=None)

    parser.add_argument("--start_idx", type=int, default=None)
    parser.add_argument("--end_idx", type=int, default=None)
    parser.add_argument("--dataset", type=str, default=config.DATASET)

    return parser.parse_args()


if __name__ == '__main__':
    train(get_args())
