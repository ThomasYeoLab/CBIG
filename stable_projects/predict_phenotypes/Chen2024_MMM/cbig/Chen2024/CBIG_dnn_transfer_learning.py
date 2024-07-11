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
import torch.nn as nn
import torch.optim as optim
from skorch import NeuralNetRegressor
from scipy.stats.stats import pearsonr
from skorch.callbacks import EarlyStopping
from sklearn.model_selection import GridSearchCV
from config import config
from CBIG_misc import cod_znormed, demean_normalize, misc_z_norm, split_tra_val_with_y, read_phes


def modify_net(net, args):
    '''
    Modify the network structure for transfer learning.

    Args:
        net (nn.Module): The network to be updated.
        args: args from command line

    Returns:
        nn.Module: Modified network.
    '''
    fc_last = net.layers[-1]
    # set number of output node (for target dataset) to be 1
    fc_last[1] = nn.Linear(fc_last[1].in_features, 1)
    # re-initializes last layer
    nn.init.xavier_uniform_(fc_last[1].weight)
    fc_last[1].bias.data.fill_(0.01)

    n_layer_for_finetune = args.n_layer_for_finetune # Number of layers to modify from the end
    start_layer_index = max(len(net.layers) - n_layer_for_finetune, 0)

    # Set requires_grad for the selected layers
    for i, layer in enumerate(net.layers):
        if i >= start_layer_index:
            layer.requires_grad = True
        else:
            layer.requires_grad = False

    return net


def finetune(x_test_k, x_test_remain, y_test_k, opt_index, device, args):
    '''finetune the DNN

    Args:
        x_test_tra (ndarray): training test FC data
        x_test_val (ndarray): validation test FC data
        x_test_remain (ndarray): remaining test FC data
        y_test_tra (ndarray): training test phenotype data
        y_test_val (ndarray): validation test phenotype data
        opt_index (int): optimal epoch for base model trained on training meta-set
        device (torch.device): the gpu that torch runs on
        args: args from command line

    Returns:
        ndarray: finetuned result, if result is worse than original, return
            none.
    '''
    k = x_test_k.shape[0]
    batch_size = k if k < 32 else 32

    y_test_k = torch.from_numpy(y_test_k.astype(np.float32).reshape(-1, 1))
    x_test_k = torch.from_numpy(x_test_k.astype(np.float32))
    x_test_remain = torch.from_numpy(x_test_remain.astype(np.float32))

    dir_model = os.path.join(args.model_dir, 'trained_model_' + args.src_dataset)
    dir_model_temp = os.path.join(dir_model, 'dnn_model_save_base')
    model_path = os.path.join(dir_model_temp, 'CBIG_dnn_epoch_' + str(opt_index) + '.pkl_torch')
    net = torch.load(model_path)
    net = modify_net(net, args)
    # early stopping by validation set (20% * K)
    early_stopping = EarlyStopping(monitor='valid_loss', patience=5, threshold_mode='rel', lower_is_better=True)
    NNR = NeuralNetRegressor(
        module=net,
        lr=1e-3,
        batch_size=batch_size,
        criterion=nn.MSELoss(),
        optimizer=optim.SGD,
        optimizer__momentum=args.momentum,
        optimizer__weight_decay=args.weight_decay,
        iterator_train__shuffle=True,
        max_epochs=args.epochs,
        callbacks= [early_stopping],
        device=device)
    # use cross-validation to tune learning rate
    params = {'lr': [1e-3, 1e-4, 1e-5]}
    gs = GridSearchCV(NNR, params, cv=5, refit=True)
    gs.fit(x_test_k, y_test_k)
    y_test_remain_pred = gs.predict(x_test_remain)
    y_test_k_pred = gs.predict(x_test_k)

    return np.squeeze(y_test_remain_pred), np.squeeze(y_test_k_pred)


def transfer_learning_dnn(y_test, x_test, opt_index, device, split_ind, args):
    '''transfer learning

    Args:
        y_test (ndarray): original test data. Shape: (#subjects, )
        x_test (ndarray): original test FC matrices. Shape: (#subjects, #roi * (#roi + 1) / 2)
        opt_index (int): optimal epoch for base model trained on training meta-set
        device (torch.device): the gpu that torch runs on
        split_ind (ndarray): array that indicate K subjects
        args: args from command line

    Returns:
        Tuple: result (normal DNN and finetuned DNN correlation and COD) of
            transfer learning
    '''
    # get split from krr pure baseline
    split = np.squeeze(split_ind)
    split_k = split == False
    split_tes = split == True
    split_tra, split_val = split_tra_val_with_y(split, y_test)

    assert np.array_equal(split_k, split_tra + split_val)
    x_test_tra = x_test[split_k]
    x_test_remain = x_test[split_tes, :]

    # z norm based on training set (K-shot)
    _, y_test, _, _ = misc_z_norm(y_test[split_k], y_test)
    y_test_tra = y_test[split_k]
    y_test_remain = y_test[split_tes]

    # exclude nan value from y_test_remain
    real_index = ~np.isnan(y_test_remain)
    y_test_remain = y_test_remain[real_index]
    x_test_remain = x_test_remain[real_index, :]

    y_pred_tuned, y_pred_k = finetune(x_test_tra, x_test_remain, y_test_tra, opt_index, device, args)
    res_cor = pearsonr(y_test_remain, y_pred_tuned)[0]
    res_cod = cod_znormed(y_test_remain, y_pred_tuned)

    return res_cor, res_cod, y_pred_k, split_k


def main(args):
    '''main function for transfer learning

    Args:
        args: args from command line

    Returns:
        None
    '''

    print('\nCBIG transfer learning with argument: ' + str(args))

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    ks = config.KS
    n_rng = config.N_RNG

    if args.exp_dataset:
        args.in_dir = config.IN_DIR_UT if args.unit_test else config.IN_DIR_EXP
        args.out_dir = config.OUT_DIR_UT if args.unit_test else config.OUT_DIR_EXP
        args.inter_dir = config.INTER_DIR_UT if args.unit_test else config.INTER_DIR_EXP
        args.model_dir = config.MODEL_DIR_UT if args.unit_test else config.MODEL_DIR_EXP
        args.src_dataset = config.DATASET_NAME_EXP['extra-large']
        args.tar_dataset = config.DATASET_NAME_EXP['test'][0]
        n_rng = config.N_RNG_EXP

    log_str = args.log_stem + '_2' + args.tar_dataset + '_result'
    meta_cor_npz = os.path.join(args.out_dir, log_str + '.npz')

    npz = np.load(os.path.join(args.inter_dir, 'dnn_base.npz'))
    tuned_by = 'cod'
    val_record = npz['val_' + tuned_by + '_record']
    temp = np.mean(val_record[:, :], axis=1)
    temp = np.convolve(temp, np.ones(3, dtype=int), 'valid') / 3
    index = np.nanargmax(temp)
    index = index + 1
    print('\nBest validation at index: ', index)
    print(np.mean(val_record[index, :]))

    npz = np.load(os.path.join(args.in_dir, args.tar_dataset,
                               args.tar_dataset + '_dnn_input.npz'), allow_pickle=True)
    y_test = npz['y_raw']
    y_test = np.array(y_test, dtype=np.float64)
    x_test = npz['x_raw']
    x_test[np.isnan(x_test)] = 0
    x_test = demean_normalize(x_test)
    x_test = np.array(x_test, dtype=np.float64)

    tes_phe = read_phes(args.in_dir, args.tar_dataset)

    # load data split
    npz = np.load(os.path.join(args.inter_dir, args.tar_dataset +
                               '_split_ind_' + str(n_rng) + '.npz'), allow_pickle=True)
    split_ind_dict = npz['split_ind_dict'].item()

    # perform meta matching with DNN and transfer learning
    start_time = time.time()
    meta_cor = np.zeros((n_rng, len(ks), len(tes_phe)))
    meta_cod = np.zeros((n_rng, len(ks), len(tes_phe)))
    if args.haufe_save:
        meta_pred_k_100 = np.zeros((n_rng, len(tes_phe), 100))
        meta_x_k_100 = np.zeros((n_rng, len(tes_phe), 100, x_test.shape[1]))
    for i in range(n_rng):
        for ik, k in enumerate(ks):
            for ib, phe in enumerate(tes_phe):
                split_ind = split_ind_dict.get(phe)[i, ik, :]
                meta_cor[i, ik, ib], meta_cod[i, ik, ib], tmp_pred, split_k = \
                    transfer_learning_dnn(y_test[:, ib], x_test, index, device, split_ind, args)
                if args.haufe_save and k == 100:
                    meta_pred_k_100[i, ib, :] = tmp_pred
                    meta_x_k_100[i, ib, :, :] = x_test[split_k, :]
        print("rng %d at %ss: cor %.5f, cod %.5f" %
              (i, time.time() - start_time, np.nanmean(meta_cor[:i + 1, :, :]),
               np.nanmean(meta_cod[:i + 1, :, :])))
        mean_cor = np.squeeze(
            np.nanmean(np.nanmean(meta_cor[:i + 1, :, :], axis=2), axis=0))
        mean_cod = np.squeeze(
            np.nanmean(np.nanmean(meta_cod[:i + 1, :, :], axis=2), axis=0))
        print(' '.join('%.6f' % tmp for tmp in mean_cor), ' COD ', ' '.join('%.6f' % tmp for tmp in mean_cod))

    np.savez(meta_cor_npz, meta_cor=meta_cor, meta_cod=meta_cod, tes_phe=tes_phe)

    if args.haufe_save:
        haufe_npz = os.path.join(args.large_data_dir,
                                 'haufe_y_pred_100_' + args.log_stem + '_2' + args.tar_dataset + '.npz')
        np.savez(haufe_npz, meta_pred_k_100=meta_pred_k_100, meta_x_k_100=meta_x_k_100)
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
    parser.add_argument('--model_dir', type=str, default=config.MODEL_DIR)
    parser.add_argument('--large_data_dir', type=str, default=config.LARGE_DATA_DIR)
    parser.add_argument('--log_stem', type=str, default='transfer_learning')
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--src_dataset', type=str, default='UKBB')
    parser.add_argument('--tar_dataset', type=str, default='HCP')

    # hyperparameter
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--momentum', type=float, default=0.9)
    parser.add_argument('--weight_decay', type=float, default=1e-5)
    parser.add_argument('--val_interval', type=int, default=1)
    parser.add_argument('--epochs', type=int, default=10)
    parser.add_argument('--n_layer_for_finetune', type=float, default=2)

    haufe_save_parser = parser.add_mutually_exclusive_group(required=False)
    haufe_save_parser.add_argument('--haufe-save', dest='haufe_save', action='store_true')
    haufe_save_parser.add_argument('--not-haufe-save', dest='haufe_save', action='store_false')
    parser.set_defaults(haufe_save=False)

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
