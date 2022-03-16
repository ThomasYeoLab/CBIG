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
from scipy.stats.stats import pearsonr

import torch
import torch.utils.data
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader

from config import config
from CBIG_model_pytorch import ukbb_multi_task_dataset
from CBIG_mics import mics_z_norm, cod_znormed, split_tra_val_with_y


def update_net(net, param, args):
    '''update the DNN for finetuning

    Args:
        net (nn.Module): the DNN to be update
        param (ndarray): index of matched phenotype
        args: args from command line

    Returns:
        nn.Module: updated net
    '''
    if args.exp_dataset:
        tmp_x = torch.from_numpy(np.random.rand(4, 190)).float().cuda()
    else:
        tmp_x = torch.from_numpy(np.random.rand(4, 1485)).float().cuda()
    # Generate a random input with shape of 4, 1485 to network
    # where 1485 is the number of element of flat 55*55 FC, 55 * 54 / 2
    tmp_y = net(tmp_x)
    if hasattr(net, 'fc4'):
        fc_last = net.fc4
        fcneg1 = 'fc4.'
        fcneg2 = 'fc3.'
    elif hasattr(net, 'fc3'):
        fc_last = net.fc3
        fcneg1 = 'fc3.'
        fcneg2 = 'fc2.'
    elif hasattr(net, 'fc2'):
        fc_last = net.fc2
        fcneg1 = 'fc2.'
        fcneg2 = 'fc1.'
    tmp_weight = fc_last[1].weight[param, :].unsqueeze(0)
    tmp_bias = fc_last[1].bias[param].unsqueeze(0)
    fc_last[1] = nn.Linear(fc_last[1].in_features, 1)
    fc_last[1].weight = nn.Parameter(tmp_weight)
    fc_last[1].bias = nn.Parameter(tmp_bias)

    tmp_y2 = net(tmp_x)
    if not np.allclose(tmp_y[:, param].data.cpu().numpy(),
                       tmp_y2.data.cpu().numpy().squeeze()):
        raise Exception('errors during update net')

    if args.n_tuned_layer == 1:
        for name, i in net.named_parameters():
            if name.startswith(fcneg1):
                i.requires_grad = True
            else:
                i.requires_grad = False
    elif args.n_tuned_layer == 2:
        for name, i in net.named_parameters():
            if name.startswith(fcneg1) or name.startswith(fcneg2):
                i.requires_grad = True
            else:
                i.requires_grad = False
    else:
        raise Exception('to many tuned layer')

    return net


def fine_tune(x_test_tra, x_test_val, x_test_remain, y_test_tra, y_test_val,
              y_test_remain, opt_index, best_phe, device, args):
    '''finetune the DNN

    Args:
        x_test_tra (ndarray): training test FC data
        x_test_val (ndarray): validation test FC data
        x_test_remain (ndarray): remaining test FC data
        y_test_tra (ndarray): training test phenotype data
        y_test_val (ndarray): validation test phenotype data
        y_test_remain (ndarray): remaining test phenotype data
        opt_index (int): optimal epoch for base model trained on training
            meta-set
        best_phe (int): index of matched phenotype
        device (torch.device): the gpu that torch runs on
        args: args from command line

    Returns:
        ndarray: finetuned result, if result is worse than original, return
            none.
    '''

    k = x_test_tra.shape[0]
    batch_size = k if k < 32 else 32

    dset_train = ukbb_multi_task_dataset(x_test_tra, y_test_tra, True)
    trainloader = DataLoader(
        dset_train, batch_size=batch_size, shuffle=True, num_workers=8)
    dset_valid = ukbb_multi_task_dataset(x_test_val, y_test_val, True)
    validLoader = DataLoader(
        dset_valid, batch_size=batch_size, shuffle=True, num_workers=8)
    dset_test = ukbb_multi_task_dataset(x_test_remain, y_test_remain, True)
    testLoader = DataLoader(
        dset_test, batch_size=batch_size, shuffle=False, num_workers=8)
    if args.exp_dataset:
        weight_path = os.path.join(
            args.out_dir, 'trained_model_ukbb',
            'dnn_model_save_base_exp_dataset',
            'CBIG_ukbb_dnn_run_0_epoch_' + str(opt_index) + '.pkl_torch')
    else:
        weight_path = os.path.join(
            args.out_dir, 'trained_model_ukbb', 'dnn_model_save_base',
            'CBIG_ukbb_dnn_run_0_epoch_' + str(opt_index) + '.pkl_torch')
    net = torch.load(weight_path)
    tmp_path = os.path.join(
        args.out_dir,
        'tmp/transfer_tmp_model_save_gpu_' + str(args.gpu) + '.pkl_torch')

    criterion = nn.MSELoss()
    optimizer = optim.SGD(
        net.parameters(),
        lr=args.lr,
        momentum=args.momentum,
        weight_decay=args.weight_decay)
    net = update_net(net, best_phe, args)

    net.train(False)
    record_loss = 0.0
    for (x, y) in validLoader:
        x, y = x.to(device), y.to(device)
        outputs = net(x)
        loss = criterion(outputs, y)
        record_loss += loss.item()
    torch.save(net.state_dict(), tmp_path)
    valid_loss = record_loss / len(validLoader)
    is_finetune_good = False

    for epoch in range(args.epochs):
        net.train(True)
        for (x, y) in trainloader:
            x, y = x.to(device), y.to(device)
            optimizer.zero_grad()
            outputs = net(x)
            loss = criterion(outputs, y)
            loss.backward()
            optimizer.step()

        if epoch % args.val_interval == args.val_interval - 1:
            net.train(False)
            record_loss = 0.0
            for (x, y) in validLoader:
                x, y = x.to(device), y.to(device)
                outputs = net(x)
                loss = criterion(outputs, y)
                record_loss += loss.item()
            record_loss = record_loss / len(validLoader)

            if record_loss < valid_loss:
                is_finetune_good = True
                valid_loss = record_loss
                torch.save(net.state_dict(), tmp_path)

    if is_finetune_good:
        net.load_state_dict(torch.load(tmp_path))
        net.train(False)
        record_pred = np.zeros((0, 1))
        for (x, y) in testLoader:
            x, y = x.to(device), y.to(device)
            outputs = net(x)
            record_pred = np.concatenate(
                (record_pred, outputs.data.cpu().numpy()), axis=0)
        return np.squeeze(record_pred)
    else:
        # print("not good")
        return None


def mm_dnn_finetune(y_test, y_pred, x_test, opt_index, device, split_ind,
                    args):
    '''meta-matching with DNN finetune

    Args:
        y_test (ndarray): original test data
        y_pred (ndarray): predicted for test subjects with base model trained
            on training meta-set
        x_test (ndarray): original test FC data
        opt_index (int): optimal epoch for base model trained on training
            meta-set
        device (torch.device): the gpu that torch runs on
        split_ind (ndarray): array that indicate K subjects
        args: args from command line

    Returns:
        Tuple: result (normal DNN and finetuned DNN correlation and COD) of
            meta-matching with DNN and finetune, and best matched phenotypes
    '''

    # get split from krr pure baseline
    split = np.squeeze(split_ind)
    split_k = split == 0
    split_tes = split == 1
    split_tra, split_val = split_tra_val_with_y(split, y_test)

    assert np.array_equal(split_k, split_tra + split_val)
    x_test_tra = x_test[split_tra, :]
    x_test_val = x_test[split_val, :]
    x_test_remain = x_test[split_tes, :]
    y_pred_k = y_pred[split_k, :]
    y_pred_remain = y_pred[split_tes, :]

    # z norm based on y_train
    _, y_test, _, _ = mics_z_norm(y_test[split_k], y_test)
    y_test_k = y_test[split_k]
    y_test_tra = y_test[split_tra]
    y_test_val = y_test[split_val]
    y_test_remain = y_test[split_tes]

    # exclude nan value from y_test_remain
    real_index = ~np.isnan(y_test_remain)
    y_test_remain = y_test_remain[real_index]
    y_pred_remain = y_pred_remain[real_index, :]
    x_test_remain = x_test_remain[real_index, :]

    best_score = float("-inf")
    best_phe = -1
    best_sign = 1
    for i in range(y_pred.shape[1]):
        sign = 1
        tmp = cod_znormed(y_test_k, y_pred_k[:, i])
        y_pred_k_tmp = -y_pred_k[:, i]
        tmp1 = cod_znormed(y_test_k, y_pred_k_tmp)
        if tmp1 > tmp:
            tmp = tmp1
            sign = -1
        if tmp >= best_score:
            best_score = tmp
            best_phe = i
            best_sign = sign
    y_mm_pred = best_sign * y_pred_remain[:, best_phe]
    res_cor = pearsonr(y_test_remain, y_mm_pred)[0]
    res_cod = cod_znormed(y_test_remain, y_mm_pred)

    y_pred_tuned = fine_tune(x_test_tra, x_test_val, x_test_remain,
                             best_sign * y_test_tra, best_sign * y_test_val,
                             best_sign * y_test_remain, opt_index, best_phe,
                             device, args)

    if y_pred_tuned is None:
        res_cor_new = res_cor
        res_cod_new = res_cod
    else:
        y_mm_pred_tuned = best_sign * y_pred_tuned
        res_cor_new = pearsonr(y_test_remain, y_mm_pred_tuned)[0]
        res_cod_new = cod_znormed(y_test_remain, y_mm_pred_tuned)

    return res_cor, res_cor_new, res_cod, res_cod_new, best_phe


def main(args):
    '''main function for meta-matching (DNN finetune)

    Args:
        args: args from command line

    Returns:
        None
    '''

    print('\nCBIG meta-matching (DNN finetune) with argument: ' + str(args))

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
    log_str = args.log_stem + '_result_' + args.split + '_fine_tune'
    if args.exp_dataset:
        log_str += '_exp_dataset'
    meta_cor_npz = os.path.join(args.out_dir,
                                log_str + '_gpu_' + str(args.gpu) + '.npz')
    meta_temp_npz = os.path.join(
        args.out_dir, 'tmp/' + log_str + '_rng_gpu_' + str(args.gpu) + '.npz')
    os.makedirs(os.path.join(args.out_dir, 'tmp'), exist_ok=True)
    ks = [10, 20, 50, 100, 200]
    n_rng = args.rng

    # load base prediction result
    if args.exp_dataset:
        npz = os.path.join(args.out_dir, 'dnn_exp_dataset_base.npz')
    else:
        npz = os.path.join(args.out_dir, 'dnn_base.npz')
    npz = np.load(npz)
    tes_res_record = npz['tes_res_record']
    val_record = npz['val_' + tuned_by + '_record']
    temp = np.mean(val_record[0, :, :], axis=1)
    temp = np.convolve(temp, np.ones(3, dtype=int), 'valid') / 3
    index = np.nanargmax(temp)
    index = index + 1
    print('\nBest validation at index: ', index)
    y_pred = tes_res_record[0, index, :, :]

    # load original data
    if args.exp_dataset:
        npz = os.path.join(args.in_dir, 'exp_dnn_input_test.npz')
    else:
        npz = os.path.join(args.in_dir, 'ukbb_dnn_input_test.npz')
    npz = np.load(npz, allow_pickle=True)
    x_test = npz['x_test']
    y_test = npz['y_test_different_set_phe']
    tes_phe = npz['tes_phe']
    tra_phe = npz['tra_phe']

    # load data split
    if args.exp_dataset:
        npz = np.load(
            os.path.join(args.inter_dir,
                         'exp_split_ind_' + str(n_rng) + '.npz'),
            allow_pickle=True)
    else:
        npz = np.load(
            os.path.join(args.inter_dir,
                         'ukbb_split_ind_' + str(n_rng) + '.npz'),
            allow_pickle=True)
    split_ind_dict = npz['split_ind_dict'].item()

    # perform meta matching with DNN and transfer learning
    start_time = time.time()
    meta_cor = np.zeros((n_rng, len(ks), len(tes_phe)))
    meta_cor_tuned = np.zeros((n_rng, len(ks), len(tes_phe)))
    meta_cod = np.zeros((n_rng, len(ks), len(tes_phe)))
    meta_cod_tuned = np.zeros((n_rng, len(ks), len(tes_phe)))
    meta_phe = np.zeros((n_rng, len(ks), len(tes_phe)))
    for i in range(n_rng):
        for ik, _ in enumerate(ks):
            for ib, phe in enumerate(tes_phe):
                split_ind = split_ind_dict.get(phe)[i, ik, :]
                meta_cor[i, ik, ib], meta_cor_tuned[i, ik, ib], meta_cod[
                    i, ik, ib], meta_cod_tuned[
                        i, ik, ib], meta_phe[i, ik, ib] = mm_dnn_finetune(
                            y_test[:, ib], y_pred, x_test, index, device,
                            split_ind, args)
        print("rng %d at %ss: cor %.5f, tuned %.5f, cod %.5f, tuned %.5f" %
              (i, time.time() - start_time, np.nanmean(meta_cor[:i + 1, :, :]),
               np.nanmean(meta_cor_tuned[:i + 1, :, :]),
               np.nanmean(meta_cod[:i + 1, :, :]),
               np.nanmean(meta_cod_tuned[:i + 1, :, :])))
        mean_cor_tuned = np.squeeze(
            np.nanmean(
                np.nanmean(meta_cor_tuned[:i + 1, :, :], axis=2), axis=0))
        mean_cod_tuned = np.squeeze(
            np.nanmean(
                np.nanmean(meta_cod_tuned[:i + 1, :, :], axis=2), axis=0))
        print(' '.join('%.6f' % tmp for tmp in mean_cor_tuned), ' COD ',
              ' '.join('%.6f' % tmp for tmp in mean_cod_tuned))
        np.savez(
            meta_temp_npz,
            meta_cor=meta_cor,
            meta_cor_tuned=meta_cor_tuned,
            meta_cod=meta_cod,
            meta_cod_tuned=meta_cod_tuned,
            current_rng=i,
            meta_phe=meta_phe,
            tes_phe=tes_phe,
            tra_phe=tra_phe)

    np.savez(
        meta_cor_npz,
        meta_cor=meta_cor,
        meta_cor_tuned=meta_cor_tuned,
        meta_cod=meta_cod,
        meta_cod_tuned=meta_cod_tuned,
        meta_phe=meta_phe,
        tes_phe=tes_phe,
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
    parser.add_argument('--inter_dir', type=str, default=config.INTER_DIR)
    parser.add_argument('--log_stem', type=str, default='meta')
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--split', type=str, default='test')
    parser.add_argument('--rng', type=int, default=100)

    exp_dataset_parser = parser.add_mutually_exclusive_group(required=False)
    exp_dataset_parser.add_argument(
        '--exp-dataset', dest='exp_dataset', action='store_true')
    exp_dataset_parser.add_argument(
        '--not-exp-dataset', dest='exp_dataset', action='store_false')
    parser.set_defaults(exp_dataset=False)

    # hyperparameter
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--momentum', type=float, default=0.9)
    parser.add_argument('--weight_decay', type=float, default=1e-7)
    parser.add_argument('--val_interval', type=int, default=10)
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--n_tuned_layer', type=int, default=2)

    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
