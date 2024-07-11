#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Chen Pansheng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import torch
import numpy as np
from scipy.stats.stats import pearsonr
from sklearn.preprocessing import normalize


def demean_normalize(val):
    '''de-mean and normalize data

    Args:
        val (ndarray): value to be de-meaned and normalized

    Returns:
        ndarray: de-meaned and normalized data
    '''
    mu = np.nanmean(val, axis=1, keepdims=True)
    val = val - mu
    val = normalize(val, axis=1, norm='l2')
    return val


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


def misc_z_norm(train_y, valid_y=None, test_y=None):
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
    if valid_y is not None:
        valid_y = valid_y - t_mu
    if test_y is not None:
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
    if valid_y is not None:
        valid_y = valid_y / t_sigma_d
    if test_y is not None:
        test_y = test_y / t_sigma_d

    # return processed y and std for future MAE calculation
    return [train_y, valid_y, test_y, t_sigma]


def read_phes(folder, dataset):
    '''read the phenotypes in a dataset

    Args:
        folder (str): folder path that stores data input files
        dataset (str): name of the dataset

    Returns:
        List: list of phenotype name, with shape of (#phenotypes, )
    '''
    txt_path = os.path.join(folder, dataset, dataset + '_phe_list.txt')
    phes = np.genfromtxt(txt_path, dtype=str)
    return phes


def get_phe_num(folder, dataset):
    '''get the number of phenotypes in a dataset

    Args:
        folder (str): folder path that stores data input files
        dataset (str): name of the dataset

    Returns:
        Int: number of phenotypes in the dataset
    '''
    subj_list_path = os.path.join(folder, dataset, dataset + '_phe_list.txt')
    with open(subj_list_path, 'r') as file:
        lines = file.readlines()
    n_phe = len(lines)
    return n_phe


def get_subj_num(folder, dataset):
    '''get the number of subjects in a dataset

    Args:
        folder (str): folder path that stores data input files
        dataset (str): name of the dataset

    Returns:
        Int: number of subjects in the dataset
    '''
    subj_list_path = os.path.join(folder, dataset, dataset + '_subj_list.txt')
    with open(subj_list_path, 'r') as file:
        lines = file.readlines()
    n_phe = len(lines)
    return n_phe


def compare_dicts(dict1, dict2, rtol=1e-05, atol=1e-08):
    """
    Compare two dictionaries element-wise with tolerance.

    Args:
        dict1 (dict): First dictionary to compare.
        dict2 (dict): Second dictionary to compare.
        rtol (float): The relative tolerance parameter (see numpy.allclose).
        atol (float): The absolute tolerance parameter (see numpy.allclose).

    Returns:
        bool: True if all elements in the dictionaries are equal within the given tolerance, False otherwise.
    """
    # Check if the dictionaries have the same keys
    if dict1.keys() != dict2.keys():
        return False

    # Iterate over the keys
    for key in dict1:
        # Check if the values are array-like and numeric
        if isinstance(dict1[key], np.ndarray) and np.issubdtype(dict1[key].dtype, np.number) and \
           isinstance(dict2[key], np.ndarray) and np.issubdtype(dict2[key].dtype, np.number):
            # Compare the arrays corresponding to each key with a tolerance for error
            if not np.allclose(dict1[key], dict2[key], rtol=rtol, atol=atol):
                return False
        else:
            # If the values are not array-like and numeric, compare them directly
            if not np.array_equal(dict1[key], dict2[key]):
                return False

    # All arrays are equal within the given tolerance
    return True



def misc_infer_metric(dataloader,
                      net,
                      criterion,
                      device,
                      t_sigma=None,
                      need_value=False,
                      output_size=1):
    '''performance inference with net on data from dataloader and calculate metric

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
    for (x, y) in dataloader:
        x, y = x.to(device), y.to(device)
        outputs = net(x)
        loss = criterion(outputs, y)
        record_loss += loss.item()
        record_real = np.concatenate((record_real, y.data.cpu().numpy()), axis=0)
        record_pred = np.concatenate((record_pred, outputs.data.cpu().numpy()), axis=0)
        if t_sigma is None:
            _, predicted = torch.max(outputs.data, 1)
            record_total += y.size(0)
            record_correct += (predicted == y.data).sum()

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


def misc_log(model_name, out_path, metric='cor', index=None, **kwargs):
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
        temp = np.mean(val_record, axis=1)
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
          np.nanmean(np.nanmean(val_cor_record[index, :], axis=0)),
          ', COD:', np.nanmean(
              np.nanmean(val_cod_record[index, :], axis=0)), ', MAE:',
          np.nanmean(np.nanmean(val_mae_record[index, :], axis=0)))

    return


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
