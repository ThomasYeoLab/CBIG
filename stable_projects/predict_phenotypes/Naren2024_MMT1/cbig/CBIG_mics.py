#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
import torch
from nipype.interfaces.fsl import AvScale


def nanpearsonr(real, pred):
    '''Compute Pearson's correlation, omit NAN

    Args:
        real (ndarray): original value
        pred (ndarray): predicted value

    Returns:
        ndarray: Correlation value
    '''

    n = real.shape[1]
    res = np.zeros((n))
    for i in range(n):
        tmp = np.logical_not(np.isnan(real[:, i]))
        res[i] = pearsonr(real[tmp, i], pred[tmp, i])[0]
    return res


def cod_znormed(real, pred):
    '''Compute COD (Coefficient of Determination) for z normed data

    Args:
        real (ndarray): original value
        pred (ndarray): predicted value

    Returns:
        float: COD value
    '''

    tot = np.sum((real)**2, axis=-1)
    res = np.sum((real - pred)**2, axis=-1)
    return 1 - res / tot


def nancod_znormed(real, pred):
    '''Compute COD (Coefficient of Determination) for z normed data, omit NAN

    Args:
        real (ndarray): original value
        pred (ndarray): predicted value

    Returns:
        ndarray: COD value
    '''

    tot = np.nansum((real)**2, axis=-2)
    res = np.nansum((real - pred)**2, axis=-2)
    return 1 - np.divide(res, tot)


def split_tra_val_with_y(split, y_test):
    '''split K subjects into training and validation

    Args:
        split (ndarray): array that indicate K subjects
        y_test (ndarray): test data to avoid same value after split train
        and validation data

    Returns:
        Tuple: split array for training and validation data
    '''

    n = split.shape[0]
    n_rng = 1
    while n_rng != 0:
        k = np.where(split == 0)[0]
        m = k.shape[0]
        m_tra = int(m * 0.8)
        k = np.random.permutation(k)
        split_tra = np.zeros((n))
        split_tra[k[:m_tra]] = 1
        split_tra = split_tra.astype(bool)
        split_val = np.zeros((n))
        split_val[k[m_tra:]] = 1
        split_val = split_val.astype(bool)
        y_test_tra = y_test[split_tra]
        if np.unique(y_test_tra).shape[0] > 1:
            n_rng = 0
        else:
            np.random.seed(100 + n_rng)
            n_rng += 1

    return split_tra, split_val


def mics_z_norm(train_y, valid_y, test_y=None):
    '''z normalize y of training, validation and test set based on training set

    Args:
        train_y (ndarray): training y data
        valid_y (ndarray): validation y data
        test_y (ndarray, optional): testing y data

    Returns:
        Tuple: contains z-normed y data and std of training y data
    '''

    # subtract mean of y of training set
    t_mu = np.nanmean(train_y, axis=0, keepdims=True)
    train_y = train_y - t_mu
    valid_y = valid_y - t_mu
    if test_y:
        test_y = test_y - t_mu

    # divide std of y of training set
    t_sigma = np.nanstd(train_y, axis=0)
    if train_y.ndim == 2:
        t_sigma_d = t_sigma[np.newaxis, :]
    else:
        t_sigma_d = t_sigma
        if t_sigma == 0:
            print('t_sigma is 0, pass divide std')
            return [train_y, valid_y, test_y, t_sigma]
    train_y = train_y / t_sigma_d
    valid_y = valid_y / t_sigma_d
    if test_y:
        test_y = test_y / t_sigma_d

    return [train_y, valid_y, test_y, t_sigma]


def mics_infer_metric(dataloader,
                      net,
                      criterion,
                      device,
                      t_sigma=None,
                      need_value=False,
                      output_size=1):
    '''performance inference with net on data from dataloader and calculate
       metric

    Args:
        dataloader: dataloader to load data for PyTorch framework
        net: PyTorch deep learning network
        criterion: criterion for loss calculation
        device: torch device indicate which GPU is running
        t_sigma (float, optional): std of training y data, only use if sex is
          not the behavioral measuers
        need_value (bool, optional): whether return record of real and
          predicted value
        output_size (int, optional): size of network output

    Returns:
        Tuple: if t_sigma is not None, correlation, MAE and loss are returned.
           If t_sigma is None, auccuracy and loss are returned. If need_value
          set to True, tuple returned also returns record of real and
          predicted y value alongside the metrics. If need_value is false,
           only metrics are returned.
    '''

    # initialize variable for record
    record_loss = 0.0
    if t_sigma is None:
        record_correct = 0.0  # count of correct prediction
        record_total = 0.0  # count of total prediction
        record_real = np.zeros((0))
        record_pred = np.zeros((0, 2))
    else:
        record_real = np.zeros((0, output_size))  # real value
        record_pred = np.zeros((0, output_size))  # prediction value

    # perform inference
    with torch.no_grad():
        for (x, y, icv) in dataloader:
            x, y, icv = x.to(device), y.to(device), icv.to(device)
            outputs = net(x, icv)
            loss = criterion(outputs, y)
            record_loss += loss.item()
            record_real = np.concatenate((record_real, y.data.cpu().numpy()),
                                         axis=0)
            record_pred = np.concatenate(
                (record_pred, outputs.data.cpu().numpy()
                 if len(outputs.data.cpu().numpy().shape) == 2 else
                 np.expand_dims(outputs.data.cpu().numpy(), axis=0)),
                axis=0)

            if t_sigma is None:
                _, predicted = torch.max(outputs.data, 1)
                record_total += y.size(0)
                record_correct += (predicted == y.data).sum()

            del loss, outputs, x, y, icv
            torch.cuda.empty_cache()

    # metric calculation
    loss = record_loss / len(dataloader)
    if t_sigma is None:
        aucc = record_correct.to(torch.float) / record_total
        if need_value:
            return aucc, loss, record_real, record_pred
        else:
            return aucc, loss
    else:
        corr = nanpearsonr(record_real, record_pred)
        cod = nancod_znormed(record_real, record_pred)
        mae = np.nanmean(np.abs(record_real - record_pred), 0) * t_sigma
        if need_value:
            return corr, cod, mae, loss, record_real, record_pred
        else:
            return corr, cod, mae, loss


def mics_log(model_name, out_path, metric='cor', index=None, **kwargs):
    '''function to calculate the final result and save the record

     Args:
        model_name (str): name of network/model
        out_path (str): path to save the log
        metric (str, optional): metric to select best validation
        index (int, optional): index of optimal epoch
        **kwargs: record of training, validation and test value

    Returns:
        None
    '''

    if index is None:
        val_record = kwargs['val_' + metric + '_record']
        temp = np.mean(np.mean(val_record, axis=0), axis=1)
        temp = np.convolve(temp, np.ones(3, dtype=int), 'valid') / 3
        index = np.nanargmax(temp)
        index = index + 1
        print('\nBest validation at index: ', index)

    val_cor_record = kwargs['val_cor_record']
    val_cod_record = kwargs['val_cod_record']
    val_mae_record = kwargs['val_mae_record']

    # save record value for future use
    file_str = model_name + '_base.npz'
    name_str = os.path.join(out_path, file_str)
    os.makedirs(out_path, exist_ok=True)
    np.savez(name_str, **kwargs)
    print('file saved:', name_str)

    # get average result for validation and test data
    print('Average validation corr:',
          np.nanmean(np.nanmean(val_cor_record[:, index, :], axis=0)),
          ', COD:', np.nanmean(np.nanmean(val_cod_record[:, index, :],
                                          axis=0)), ', MAE:',
          np.nanmean(np.nanmean(val_mae_record[:, index, :], axis=0)))

    return


def save_model(net, file_trained, pre_index, metric='cor', **kwargs):
    '''function to save model

    Args:
        net : PyTorch deep learning network
        file_trained (str): path to save the net
        pre_index (int): index of optimal epoch
        metric (str): metric to select best validation
        **kwargs: record of training, validation and test value

    Returns:
        pre_index (int): index of optimal epoch
    '''

    val_record = kwargs['val_' + metric + '_record']
    temp = np.mean(np.mean(val_record, axis=0), axis=1)
    temp = np.convolve(temp, np.ones(3, dtype=int), 'valid') / 3
    index = np.nanargmax(temp)
    index = index + 1
    if pre_index != index:
        pre_index = index
        torch.save(net, file_trained)

    return pre_index


def print_result(tra_cor, tra_cod, tra_mae, val_los, val_cor, val_cod, val_mae,
                 epoch):
    '''print the result during training

    Args:
        tra_cor (ndarray): correlation of training data
        tra_cod (ndarray): COD of training data
        tra_mae (ndarray): MAE of training data
        val_los (ndarray): loss of validation data
        val_cor (ndarray): correlation of validation data
        val_cod (ndarray): COD of validation data
        val_mae (ndarray): MAE of validation data
        epoch (int): epoch number

    Returns:
        None
    '''

    print('Epoch', epoch, 'train cor %.6f' % tra_cor.mean(),
          'cod %.6f' % tra_cod.mean(), 'mae %.6f' % tra_mae.mean(),
          'valid cor %.6f' % val_cor.mean(), 'cod %.6f' % val_cod.mean(),
          'mae %.6f' % val_mae.mean(), 'loss %.6f' % val_los)


def read_datapath(data_dir, sub_list, dataset):
    '''read data path for dataloader

    Args:
        data_dir (str): directory of data
        sub_list (list): subject list of each dataset
        dataset (str): name of dataset
    Returns:
        data (ndarray): data path of all subjects in a dataset
    '''

    data = []
    if dataset == 'ukbb':
        for tmp in sub_list:
            data.append(
                os.path.join(data_dir, tmp + '_T1_MNI1mm_linear.nii.gz'))
    elif dataset == 'HCPYA' or dataset == 'HCPA' or 'example_test_dataset' \
        or 'unit_tests_test_dataset' or 'example_train_dataset' \
            or 'unit_tests_train_dataset':
        for tmp in sub_list:
            data.append(
                os.path.join(data_dir, tmp, 'T1_MNILinear_1mm_newsize.nii.gz'))

    return np.array(data)


def load_icv(icv, tra_sub_list, args):
    '''load icv for SFCN model

    Args:
        icv (ndarray): icv data
        tra_sub_list (list): subject list for training
        args (argparse.ArgumentParser) : args that could be used by
          other function

    Returns:
        tra_icv (ndarray): icv data for training subjects
    '''

    tra_mask_icv = icv['eid'].isin(tra_sub_list[0].values.tolist())
    tra_icv_raw = icv.loc[tra_mask_icv, ['inverse_determinant']]
    tra_icv_raw = tra_icv_raw.values.astype(float)
    train_val_split_ukbb = np.load(os.path.join(
        args.model_dir, 'training_process/train_val_split_ukbb.npz'),
                                   allow_pickle=True)
    split_tra_ukbb = train_val_split_ukbb['split_tra']
    tra_icv = tra_icv_raw[split_tra_ukbb, :]

    return tra_icv


def load_data(args):
    '''load data used for training or finetune

    Args:
        args: args from command line

    Returns:
        icv_test (ndarray): icv data in meta-test set
        tes_sub_list (ndarray): subject list in meta-test set
        tes_phe_name (ndarray): phenotype list in meta-test set
        y_test (ndarray): output data in meta-test set

    '''

    tes_phe = pd.read_csv(os.path.join(args.phe_dir))
    icv = pd.read_csv(os.path.join(args.icv_dir))
    tes_sub_list = pd.read_table(os.path.join(args.sub_dir), header=None)
    tra_sub_list = pd.read_table(os.path.join(args.tra_sub_dir), header=None)
    tes_mask = tes_phe['eid'].isin(tes_sub_list[0].values.tolist())
    tes_phe_name = tes_phe.columns[args.start_idx:args.end_idx]
    tes_phe = tes_phe.loc[tes_mask, tes_phe_name]
    tes_mask_icv = icv['eid'].isin(tes_sub_list[0].values.tolist())
    tes_icv = icv.loc[tes_mask_icv, ['inverse_determinant']]

    if args.across_dataset:
        ukbb_icv = pd.read_csv(os.path.join(args.ukbb_icv_dir))
        tra_icv = load_icv(ukbb_icv, tra_sub_list, args)
    else:
        tra_icv = load_icv(icv, tra_sub_list, args)

    y_test = tes_phe.values.astype(float)
    tes_icv = tes_icv.values.astype(float)

    icv_train, icv_test, _, t_sigma = mics_z_norm(tra_icv, tes_icv)

    return icv_test, tes_sub_list, tes_phe_name, y_test


def extract_icv_affine(args):
    '''extract icv from acpc2MNILinear.mat used for SFCN

    Args:
        args: args from command line

    Returns:
        None
    '''

    filename = "acpc2MNILinear.mat"
    sub_list = pd.read_table(args.sub_path, header=None)
    sub_list = sub_list.values.squeeze().tolist()

    failed_subjects = []

    sub_id_list = []
    determinant = []
    inverse_determinant = []
    avscale = AvScale()

    for sub_id in sub_list:
        try:
            avscale.inputs.mat_file = os.path.join(args.data_path, str(sub_id),
                                                   filename)
            res = avscale.run()
            d = res.outputs.determinant
            sub_id_list.append(sub_id)
            determinant.append(d)
            inverse_determinant.append(round(1 / d, 4))
        except KeyError:
            failed_subjects.append(sub_id)

    df = pd.DataFrame(list(zip(sub_id_list, determinant, inverse_determinant)),
                      columns=['eid', 'determinant', 'inverse_determinant'])
    df.to_csv(args.save_path, index=None)
