#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import random
import argparse
import numpy as np
from scipy.stats.stats import pearsonr
from skorch.callbacks import EarlyStopping
from skorch import NeuralNetRegressor
from sklearn.model_selection import GridSearchCV

import torch
import torch.utils.data
import torch.nn as nn
import torch.optim as optim

from cbig.CBIG_model_pytorch import transfer_model
from cbig.CBIG_mics import mics_z_norm, cod_znormed, \
    split_tra_val_with_y, load_data
from cbig.config import config


def fine_tune(x_test_tra, x_test_remain, y_test_tra, best_phe, device,
              pretrained_net, args, icv_test_tra, icv_test_remain):
    '''finetune the DNN for meta-matching finetune

    Args:
        x_test_tra (tensor): input data in meta-test set for training
        x_test_remain (tensor): input data in meta-test set for testing
        y_test_tra (ndarray) : output data in meta-test set for training
        best_phe (ndarray) : selected best output node used for
          initializing last layer of model
        device (str): cuda or cpu
        pretrained_net : pretrained model
        args (argparse.ArgumentParser) : args that could be used by
          other functions
        icv_test_tra (ndarray): icv data in meta-test set for training
        icv_test_remain (list): icv data in meta-test set for testing

    Returns:
        record_pred (ndarray) : prediction for testing data in meta-test set
        y_test_tra_pred (ndarray) : prediction for training data in
          meta-test set

    '''

    batch_size = args.batch_size

    x_test_tra = x_test_tra.numpy()
    x_test_remain = x_test_remain.numpy()
    y_test_tra = y_test_tra.astype(np.float32).reshape(-1, 1)

    net = NeuralNetRegressor(module=transfer_model(pretrained_net,
                                                   initial_last_layer=True,
                                                   param=best_phe,
                                                   in_c=args.in_c,
                                                   out_c=args.out_c),
                             lr=args.lr,
                             batch_size=batch_size,
                             criterion=nn.MSELoss(),
                             optimizer=optim.SGD,
                             optimizer__momentum=args.momentum,
                             optimizer__weight_decay=args.weight_decay,
                             train_split=None,
                             iterator_train__shuffle=False,
                             max_epochs=args.epochs,
                             device=device)

    params = {'lr': [1e-7, 5e-7, 1e-6, 5e-6]}
    clf = GridSearchCV(net, params, cv=5, refit=False)
    inputs = np.concatenate((x_test_tra, icv_test_tra), axis=1)
    clf.fit(inputs, y_test_tra)

    opt_lr = clf.best_params_['lr']
    del clf
    torch.cuda.empty_cache()

    early_stopping = EarlyStopping(monitor='valid_loss',
                                   patience=5,
                                   threshold=0.0001,
                                   threshold_mode='rel',
                                   lower_is_better=True)
    net = NeuralNetRegressor(module=transfer_model(pretrained_net,
                                                   initial_last_layer=True,
                                                   param=best_phe,
                                                   in_c=args.in_c,
                                                   out_c=args.out_c),
                             lr=opt_lr,
                             batch_size=batch_size,
                             criterion=nn.MSELoss(),
                             optimizer=optim.SGD,
                             optimizer__momentum=args.momentum,
                             optimizer__weight_decay=args.weight_decay,
                             iterator_train__shuffle=True,
                             max_epochs=args.epochs,
                             callbacks=[early_stopping],
                             device=device)
    inputs = np.concatenate((x_test_tra, icv_test_tra), axis=1)
    net.fit(inputs, y_test_tra)
    y_test_tra_pred = net.predict(inputs)

    record_pred = np.zeros((0, 1))
    inputs = np.concatenate((x_test_remain, icv_test_remain), axis=1)
    outputs = net.predict(inputs)
    record_pred = np.concatenate((record_pred, outputs), axis=0)

    del outputs
    torch.cuda.empty_cache()

    return np.squeeze(record_pred), np.squeeze(y_test_tra_pred)


def mm_finetune_wrapper(y_test, y_pred, x_test, device, split_ind,
                        pretrained_net, args, icv):
    '''finetune wrapper for meta-matching finetune

    Args:
        y_test (ndarray): output data in meta-test set
        y_pred (ndarray): prediction in meta-test set from pretrained model
        x_test (tensor): input data in meta-test set
        device (str): cuda or cpu
        split_ind (ndarray): training and test split
        pretrained_net : pretrained model
        args (argparse.ArgumentParser) : args that could be used by
          other function
        icv (ndarray) : icv data in meta-test set

    Returns:
        res_cor (float) : correlation between prediction and ground truth
        res_cod (float) : COD between prediction and ground truth
        y_pred (ndarray) : prediction for all subjects in meta-test set

    '''

    # get split from krr pure baseline
    split = np.squeeze(split_ind)
    split_k = split == 0
    split_tes = split == 1
    split_tra, split_val = split_tra_val_with_y(split, y_test)

    assert np.array_equal(split_k, split_tra + split_val)
    x_test_tra = x_test[split_k]
    x_test_remain = x_test[split_tes]
    y_pred_k = y_pred[split_k, :]
    y_pred_remain = y_pred[split_tes, :]
    icv_test_tra = icv[split_k]
    icv_test_remain = icv[split_tes]

    # z norm based on y_train
    _, y_test, _, t_sigma = mics_z_norm(y_test[split_k], y_test)
    y_test_tra = y_test[split_k]
    y_test_remain = y_test[split_tes]

    # exclude nan value from y_test_remain
    real_index = ~np.isnan(y_test_remain)
    y_test_remain = y_test_remain[real_index]
    y_pred_remain = y_pred_remain[real_index, :]
    x_test_remain = x_test_remain[real_index]
    icv_test_remain = icv_test_remain[real_index, :]

    best_score = float("-inf")
    best_phe = -1
    best_sign = 1
    if args.metric == 'cod':
        for i in range(y_pred.shape[1]):
            sign = 1
            tmp = cod_znormed(y_test_tra, y_pred_k[:, i])
            y_pred_k_tmp = -y_pred_k[:, i]
            tmp1 = cod_znormed(y_test_tra, y_pred_k_tmp)
            if tmp1 > tmp:
                tmp = tmp1
                sign = -1
            if tmp >= best_score:
                best_score = tmp
                best_phe = i
                best_sign = sign
    else:
        for i in range(y_pred.shape[1]):
            sign = 1
            tmp = pearsonr(y_test_tra, y_pred_k[:, i])[0]
            y_pred_k_tmp = -y_pred_k[:, i]
            tmp1 = pearsonr(y_test_tra, y_pred_k_tmp)[0]
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

    y_pred_tuned, y_pred_tuned_k = fine_tune(x_test_tra, x_test_remain,
                                             best_sign * y_test_tra, best_phe,
                                             device, pretrained_net, args,
                                             icv_test_tra, icv_test_remain)

    if y_pred_tuned is None:
        res_cor_new = res_cor
        res_cod_new = res_cod
    else:
        y_mm_pred_tuned = best_sign * y_pred_tuned
        res_cor_new = pearsonr(y_test_remain, y_mm_pred_tuned)[0]
        res_cod_new = cod_znormed(y_test_remain, y_mm_pred_tuned)

    y_pred_new = np.concatenate(
        (best_sign * y_pred_tuned_k, best_sign * y_pred_tuned), axis=0)

    return res_cor_new, res_cod_new, y_pred_new


def main(args):
    '''main function for meta-matching finetune

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
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.enabled = True

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    ks = args.ks
    n_rng = args.n_rng
    dataset = args.dataset

    # load data
    icv_test, tes_sub_list, tes_phe_name, y_test = load_data(args)
    x_test = torch.load(
        os.path.join(args.data_dir, 'dnn_test_penultimate_output.pt'))

    n_expand_dim = int(x_test.shape[2] * x_test.shape[3] * x_test.shape[4])
    tmp_icv_test = np.tile(icv_test, (1, n_expand_dim))
    icv_test = tmp_icv_test.reshape((tmp_icv_test.shape[0], 1, x_test.shape[2],
                                     x_test.shape[3], x_test.shape[4]))

    log_str = args.log_stem + '_result'
    meta_cor_npz = os.path.join(args.out_dir + args.out_subdir,
                                log_str + '_test.npz')
    os.makedirs(os.path.join(args.out_dir + args.out_subdir, 'tmp'),
                exist_ok=True)
    model_path = os.path.join(args.out_dir + args.out_subdir, 'model')
    os.makedirs(model_path, exist_ok=True)

    # load base prediction result
    npz = os.path.join(args.data_dir, 'dnn_test_base.npz')
    npz = np.load(npz)
    tes_res_record = npz['tes_res_record']
    y_pred = tes_res_record[0, 0, :, :]

    # load pretrained model
    opt_index = args.index
    weight_path = os.path.join(
        args.model_dir, 'dnn_model_save_base',
        'CBIG_ukbb_dnn_run_0_epoch_' + str(opt_index) + '.pkl_torch')
    pretrained_net = torch.load(
        weight_path)  # map_location=torch.device('cpu')

    # load data split
    if args.across_dataset:
        npz = np.load(os.path.join(
            args.inter_dir, dataset + '_split_ind_' + str(n_rng) + '.npz'),
                      allow_pickle=True)
    else:
        npz = np.load(os.path.join(args.inter_dir,
                                   'ukbb_split_ind_' + str(n_rng) + '.npz'),
                      allow_pickle=True)
    split_ind_dict = npz['split_ind_dict'].item()

    meta_cor_tuned = np.zeros((n_rng, len(ks), len(tes_phe_name)))
    meta_cod_tuned = np.zeros((n_rng, len(ks), len(tes_phe_name)))
    pred = np.zeros((n_rng, len(ks), len(tes_phe_name), x_test.shape[0]))

    for i in range(n_rng):
        for ik, k in enumerate(ks):
            for ib, phe in enumerate(tes_phe_name):
                split_ind = split_ind_dict.get(phe)[i, ik, :]

                meta_cor_tuned[i, ik, ib], meta_cod_tuned[
                    i, ik, ib], tmp_pred = mm_finetune_wrapper(
                        y_test[:, ib], y_pred, x_test, device, split_ind,
                        pretrained_net, args, icv_test)

                pred[i, ik, ib, :tmp_pred.shape[0]] = tmp_pred

    np.savez(meta_cor_npz,
             meta_cor_tuned=meta_cor_tuned,
             meta_cod_tuned=meta_cod_tuned,
             pred=pred)
    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by
          other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--data_dir', type=str, default=None)
    parser.add_argument('--dataset', type=str, default=config.DATASET)
    parser.add_argument('--model_dir', type=str, default=None)
    parser.add_argument('--phe_dir', type=str, default=None)
    parser.add_argument('--icv_dir', type=str, default=None)
    parser.add_argument('--sub_dir', type=str, default=None)
    parser.add_argument('--tra_sub_dir', type=str, default=None)
    parser.add_argument('--out_dir', '-o', type=str, default=None)
    parser.add_argument('--out_subdir', '-osub', type=str, default=None)
    parser.add_argument('--inter_dir', type=str, default=None)
    parser.add_argument('--seed', type=int, default=config.SEED)
    parser.add_argument('--batch_size', type=int, default=None)
    parser.add_argument('--epochs', type=int, default=None)
    parser.add_argument('--n_rng', type=int, default=None)
    parser.add_argument('--metric', type=str, default='cod')
    parser.add_argument('--log_stem', type=str, default='meta')
    parser.add_argument('--split', type=str, default='test')
    across_dataset_parser = parser.add_mutually_exclusive_group(required=False)
    across_dataset_parser.add_argument('--across-dataset',
                                       dest='across_dataset',
                                       action='store_true')
    across_dataset_parser.add_argument('--not-across-dataset',
                                       dest='across_dataset',
                                       action='store_false')
    parser.set_defaults(across_dataset=False)

    # hyperparameter
    parser.add_argument("--ks", type=int, nargs='+', default=None)
    parser.add_argument('--lr', type=float, default=config.LR)
    parser.add_argument('--momentum', type=float, default=config.MOMENTUM)
    parser.add_argument('--weight_decay',
                        type=float,
                        default=config.WEIGHT_DECAY)

    parser.add_argument('--index', type=int, default=None)
    parser.add_argument("--start_idx", type=int, default=None)
    parser.add_argument("--end_idx", type=int, default=None)
    parser.add_argument("--ukbb_icv_dir", type=str, default=None)
    parser.add_argument("--in_c", type=int, default=None)
    parser.add_argument("--out_c", type=int, default=None)

    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
