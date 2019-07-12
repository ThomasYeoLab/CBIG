#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import time
import numpy as np
from scipy.stats.stats import pearsonr

import torch

from utils import load_graph, preprocess_adj, normalized_laplacian
from utils import rescale_laplacian, chebyshev_polynomial


def mics_z_norm(train_y, valid_y, test_y):
    '''z normalize y of training, validation and test set based on training set

    Args:
        train_y (ndarray): training y data
        valid_y (ndarray): validation y data
        test_y (ndarray): testing y data

    Returns:
        Tuple: contains z-normed y data and std of training y data
    '''
    # subtract mean of y of training set
    t_mu = train_y.mean(axis=0, keepdims=True)
    train_y = train_y - t_mu
    valid_y = valid_y - t_mu
    test_y = test_y - t_mu

    # divide std of y of training set
    t_sigma = train_y.std(axis=0)
    train_y = train_y / t_sigma[np.newaxis, :]
    valid_y = valid_y / t_sigma[np.newaxis, :]
    test_y = test_y / t_sigma[np.newaxis, :]

    # return processed y and std for future MAE calculation
    return [train_y, valid_y, test_y, t_sigma]


def mics_z_norm_gcnn(input_y, train_mask):
    """z normalize y based on training data

    Args:
        input_y (ndarray): y data
        train_mask (ndarray): mask of training data

    Returns:
        Tuple: contains z-normed y data and std of training y data
    """
    # get mean and std of training data
    y_tra = input_y[train_mask, :]
    t_mu = y_tra.mean(axis=0)
    t_sigma = y_tra.std(axis=0)

    # perform z-norm
    input_y = input_y - t_mu[np.newaxis, :]
    input_y = input_y / t_sigma[np.newaxis, :]
    return [input_y, t_sigma]


def mics_z_norm_test(train_valid_y, test_y):
    """z normalize y test set based on training data for HCP dataset

    Args:
        train_valid_y (list): list of y data for both training and validation
        test_y (ndarray): test y data

    Returns:
        Tuple: z normed test y data, and std of training y data
    """
    base_y = np.vstack(train_valid_y)
    t_v_mu = base_y.mean(axis=0)
    test_y = test_y - t_v_mu[np.newaxis, :]
    t_v_sigma = base_y.std(axis=0)
    test_y = test_y / t_v_sigma[np.newaxis, :]
    return test_y, t_v_sigma


def mics_train_valid_split(train_valid_x,
                           train_valid_y,
                           fold=None,
                           is_bnc=False):
    """split training and validation data (HCP only)

    Args:
        train_valid_x (list): list of y data for both training and validation
        train_valid_y (list): list of x data for both training and validation
        fold (int, optional): index of fold for validation, if None, no
            validation is going to be returned
        is_bnc (bool, optional): whether function is used for brainnetcnn

    Returns:
        Tuple: if fold is None, all data in list train_valid_x and y are
            combined as training x and y. If fold is not None, the
            corresponding fold is returned as validation data, while the
            remaining folds are combined as training data.
    """
    if fold is not None:
        valid_index = fold
        valid_x = train_valid_x[valid_index]
        valid_y = train_valid_y[valid_index]
        train_valid_x = np.delete(train_valid_x, valid_index, axis=0)
        train_valid_y = np.delete(train_valid_y, valid_index, axis=0)

    tmp = list(train_valid_x[0].shape)
    tmp[0] = 0
    train_x = np.zeros(tmp)
    train_y = np.zeros((0, train_valid_y[0].shape[-1]))
    for i in range(len(train_valid_x)):
        train_x = np.concatenate((train_x, train_valid_x[i]), axis=0)
        train_y = np.concatenate((train_y, train_valid_y[i]), axis=0)

    if is_bnc:
        train_x = np.expand_dims(train_x, axis=-1)

    if fold is not None:
        if is_bnc:
            valid_x = np.expand_dims(valid_x, axis=-1)
        t_mu = train_y.mean(axis=0)
        train_y = train_y - t_mu[np.newaxis, :]
        valid_y = valid_y - t_mu[np.newaxis, :]
        t_sigma = train_y.std(axis=0)
        train_y = train_y / t_sigma[np.newaxis, :]
        valid_y = valid_y / t_sigma[np.newaxis, :]
        return [train_x, valid_x, train_y, valid_y]

    t_mu = train_y.mean(axis=0)
    train_y = train_y - t_mu[np.newaxis, :]
    t_sigma = train_y.std(axis=0)
    train_y = train_y / t_sigma[np.newaxis, :]
    return [train_x, train_y]


def mics_train_valid_mask_split(train_valid_mask, fold=None):
    """split training and validation mask for gcnn (HCP only)

    Args:
        train_valid_mask (list): list of training and validation mask
        fold (int, optional): index of fold for validation, if None, no
            validation is going to be returned

    Returns:
        Tuple: training and validation mask
    """
    # Data split
    if fold is not None:
        valid_mask = train_valid_mask[fold]
        train_list = np.delete(train_valid_mask, fold, axis=0)
    else:
        valid_mask = None
        train_list = train_valid_mask

    train_mask = np.zeros(train_valid_mask[0].shape)
    for i in range(len(train_list)):
        train_mask = np.logical_or(train_mask, train_list[i])

    return [train_mask, valid_mask]


def mics_hcp_log(model_name, out_path, **kwargs):
    """calculate the test result and save the log

    Args:
        model_name (str): name of the model
        out_path (str): path to save the log npz file
        **kwargs: record of training, validation and test value

    Returns:
        None
    """
    val_cor = kwargs['val_cor_log']
    tes_cor = kwargs['tes_cor_log']
    n_folds = tes_cor.shape[0]

    temp = np.mean(val_cor, axis=-1)
    temp = np.mean(temp, axis=1)
    index = np.argmax(temp, axis=-1)
    print('Optimal index for each fold at:', index)
    result = np.array([tes_cor[i, index[i], :] for i in range(n_folds)])
    # avg = np.mean(result, axis=0)
    # err = np.std(result, axis=0) / np.sqrt(n_folds)
    temp = np.mean(result, axis=1)
    print('Optimal result for each fold:', temp)
    avg_a = np.mean(temp, axis=0)
    # err_a = np.std(temp, axis=0) / np.sqrt(n_folds)
    print('Final test result:', avg_a)
    kwargs['metric'] = avg_a

    # save record value for future use
    date_str = time.strftime("%Y_%m_%d_%H_%M")
    os.makedirs(out_path, exist_ok=True)
    file_str = 'HCP_' + model_name + '_' + date_str + '.npz'
    name_str = os.path.join(out_path, file_str)
    np.savez(name_str, **kwargs)
    print('log saved at:', file_str)

    return


def mics_hcp_infer(model, x, y, sigma, x_train=None, y_train=None):
    """evaluate model prediction for given data (HCP only)

    Args:
        model (keras.models.Model): keras DNN model
        x (ndarray): input x data
        y (ndarray): y data
        sigma (ndarray): std of training y data
        x_train (ndarray, optional): training x data
        y_train (ndarray, optional): training y data

    Returns:
        Tuple: correlation and MAE between real and predicted y, and predicted
            y value
    """
    y_pred = model.predict(x, batch_size=48, verbose=0)
    cor = np.zeros((y.shape[-1]))
    mae = np.zeros((y.shape[-1]))
    for i in range(y.shape[-1]):
        cor[i] = pearsonr(y_pred[:, i], y[:, i])[0]
        mae[i] = np.mean(np.abs(y_pred[:, i] - y[:, i])) * sigma[i]
    if x_train is None:
        return cor, mae, y_pred
    else:
        y_pred_t = model.predict(x_train, batch_size=48, verbose=0)
        cor_train = np.zeros((y_train.shape[-1]))
        for i in range(y_train.shape[-1]):
            cor_train[i] = pearsonr(y_pred_t[:, i], y_train[:, i])[0]
        return cor, mae, y_pred, cor_train


def mics_hcp_gcnn_eval(preds, input_y, mask, sigma=None, train_mask=None):
    """evaluate model prediction for given data (HCP and gcnn only)

    Args:
        preds (ndarray): predicted y value
        input_y (ndarray): real y value
        mask (ndarray): mask on y value
        sigma (ndarray, optional): std of training y data
        train_mask (ndarray, optional): mask on training y value

    Returns:
        TYPE: correlation, loss and MAE between real and predicted y
    """
    index = np.nonzero(mask)[0]
    pred = preds[index, :]
    real = input_y[index, :]

    los = np.mean(np.mean(np.square(pred - real), axis=-1))
    cor = np.zeros((input_y.shape[-1]))
    mae = np.zeros((input_y.shape[-1]))
    for i in range(input_y.shape[-1]):
        cor[i] = pearsonr(pred[:, i], real[:, i])[0]
        if sigma is not None:
            mae[i] = np.mean(np.abs(pred[:, i] - real[:, i])) * sigma[i]

    if train_mask is None:
        return cor, los, mae
    else:
        index = np.nonzero(train_mask)[0]
        pred = preds[index, :]
        real = input_y[index, :]
        cor_train = np.zeros((input_y.shape[-1]))
        for i in range(input_y.shape[-1]):
            cor_train[i] = pearsonr(pred[:, i], real[:, i])[0]
        return cor, los, mae, cor_train


def mics_infer_metric(dataloader,
                      net,
                      criterion,
                      device,
                      t_sigma=None,
                      need_value=False):
    '''performance inference with net on data from dataloader and calculate
        metric

    Args:
        dataloader: dataloader to load data for PyTorch framework
        net: PyTorch deep learning network
        criterion: criterion for loss calculation
        t_sigma (float, optional): std of training y data, only use if sex is
            not the behavioral measuers
        need_value (bool, optional): whether return record of real and
            predicted value

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
        record_real = np.zeros((0, 1))  # real value
        record_pred = np.zeros((0, 1))  # prediction value

    # perform inference
    for (x, y) in dataloader:
        x, y = x.to(device), y.to(device)
        outputs = net(x)
        loss = criterion(outputs, y)
        record_loss += loss.item()
        record_real = np.concatenate((record_real, y.data.cpu().numpy()),
                                     axis=0)
        record_pred = np.concatenate((record_pred, outputs.data.cpu().numpy()),
                                     axis=0)
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
        corr = pearsonr(record_real, record_pred)[0]
        mae = np.mean(np.abs(record_real - record_pred)) * t_sigma
        if need_value:
            return corr, mae, loss, record_real, record_pred
        else:
            return corr, mae, loss


def mics_log(model_name, out_path, index=None, item=None, **kwargs):
    '''function to calculate the final result and save the record

    Args:
        model_name (str): name of network/model
        index (int): index of optimal epoch
        out_path (str): path to save the log
        item (float, optional): indicate which behavioral meausers is predicted
        **kwargs: record of training, validation and test value

    Returns:
        None
    '''
    date_str = time.strftime("%Y_%m_%d_%H_%M")

    if index is None:
        if item is None:
            val_record = kwargs['val_auc_record']
        else:
            val_record = kwargs['val_cor_record']
        temp = np.mean(val_record, axis=0)
        temp = np.convolve(temp, np.ones(3, dtype=int), 'valid') / 3
        index = np.nanargmax(temp)
        index = index + 1
        print('\nBest validation at index: ', index)

    if item is None:
        val_auc_record = kwargs['val_auc_record']
        tes_auc_record = kwargs['tes_auc_record']
        tes_res_record = kwargs['tes_res_record']
        final_original = kwargs['final_original']

        # get result at that epoch for both validation and test
        print('Average validation aucc:',
              np.nanmean(val_auc_record[:, index], axis=0))
        print('Average test aucc:', np.nanmean(
            tes_auc_record[:, index], axis=0))

        # get ensamble result for test data
        final_predict = np.argmax(
            np.nanmean(tes_res_record[:, index, :, :], axis=0), axis=1)
        if len(final_original.shape) == 2:
            final_original = np.argmax(final_original, axis=1)
        n_test = float(final_original.shape[0])
        metric = (final_predict == final_original).sum() / n_test
        print('Final averaged test aucc', metric)

        file_str = model_name + '_sex_' + date_str + '.npz'
    else:
        val_cor_record = kwargs['val_cor_record']
        val_mae_record = kwargs['val_mae_record']
        tes_cor_record = kwargs['tes_cor_record']
        tes_mae_record = kwargs['tes_mae_record']
        tes_res_record = kwargs['tes_res_record']
        final_original = kwargs['final_original']
        t_sigma = kwargs['t_sigma']

        # get average result for validation and test data
        print('Average validation corr:',
              np.nanmean(val_cor_record[:, index], axis=0), ', MAE:',
              np.nanmean(val_mae_record[:, index], axis=0))
        print('Average test corr', np.nanmean(
            tes_cor_record[:, index], axis=0), ', MAE',
              np.nanmean(tes_mae_record[:, index], axis=0))

        # get ensamble result for test data
        final_predict = np.nanmean(tes_res_record[:, index, :], axis=0)
        final_original = np.squeeze(final_original)
        metric = pearsonr(final_predict, final_original)[0]
        print('Final ensemble test corr', metric, ', MAE',
              np.nanmean(np.abs(final_predict - final_original)) * t_sigma)

        file_str = model_name + '_pred_' + str(item) + '_' + date_str + '.npz'

    kwargs['final_predict'] = final_predict
    kwargs['metric'] = metric

    # save record value for future use
    name_str = os.path.join(out_path, file_str)
    os.makedirs(out_path, exist_ok=True)
    np.savez(name_str, **kwargs)
    print('file saved:', name_str)

    return


def mics_graph_matrix(num_subject, graph_folder, GRAPH_ADJ, FILTER,
                      MAX_DEGREE):
    """Generate graph matrix for GCNN

    Args:
        num_subject (int): number of subject for data
        graph_folder (str): location of folder for graph
        GRAPH_ADJ (str): the filename of graph
        FILTER (str): type of gcnn filter
        MAX_DEGREE (int): degree of Chebyshev polynomial

    Returns:
        Tuple: contains the graph_matrix and number of support used for GCNN

    Raises:
        Exception: invalid FILTER type
    """
    SYM_NORM = True  # symmetric (True) vs. left-only (False) normalization

    # build the graph
    A = load_graph(dimension=num_subject, path=graph_folder, graph=GRAPH_ADJ)

    # estimate the laplacian
    if FILTER == 'localpool':
        """ Local pooling filters
        (see 'renormalization trick' in Kipf & Welling, arXiv 2016)
        """
        print('Using local pooling filters...')
        A_ = preprocess_adj(A, SYM_NORM)
        support = 1
        graph_matrix = [A_]
    elif FILTER == 'chebyshev':
        """ Chebyshev polynomial basis filters
        (Defferard et al., NIPS 2016)
        """
        print('Using Chebyshev polynomial basis filters...')
        L = normalized_laplacian(A, SYM_NORM)
        L_scaled = rescale_laplacian(L)
        T_k = chebyshev_polynomial(L_scaled, MAX_DEGREE)
        support = MAX_DEGREE + 1
        graph_matrix = T_k
    else:
        raise Exception('Invalid filter type.')

    return graph_matrix, support


def mics_eval(preds, input_y, train_mask, valid_mask, test_mask, t_sigma=None):
    """evaluate the prediction for GCNN

    Args:
        preds (ndarray): GCNN prediction
        input_y (ndarray): original y data
        train_mask (ndarray): mask of training subjects
        valid_mask (ndarray): mask of validation subjects
        test_mask (ndarray): mask of testing subjects
        t_sigma (float, optional): std of training y data, only use if sex is
            not the behavioral measuers

    Returns:
        Tuple: if t_sigma is None, return the accuracy of GCNN predition. If
            t_sigma is not None, return the loss, correlation and MAE result.
    """
    val_index = np.nonzero(valid_mask)[0]
    tes_index = np.nonzero(test_mask)[0]

    if t_sigma is None:
        val_pred = np.argmax(preds[val_index, :], axis=1)
        tes_pred = np.argmax(preds[tes_index, :], axis=1)
        tra_pred = np.argmax(preds[train_mask, :], axis=1)
        val_real = np.argmax(input_y[val_index, :], axis=1)
        tes_real = np.argmax(input_y[tes_index, :], axis=1)
        tra_real = np.argmax(input_y[train_mask, :], axis=1)

        val_auc = (val_pred == val_real).mean()
        tes_auc = (tes_pred == tes_real).mean()
        tra_auc = (tra_pred == tra_real).mean()

        return [val_auc, tes_auc, tra_auc]

    val_pred = np.squeeze(preds[val_index])
    tes_pred = np.squeeze(preds[tes_index])
    tra_pred = np.squeeze(preds[train_mask])
    val_real = np.squeeze(input_y[val_index])
    tes_real = np.squeeze(input_y[tes_index])
    tra_real = np.squeeze(input_y[train_mask])

    val_los = np.mean(np.square(val_pred - val_real), axis=-1)
    tes_los = np.mean(np.square(tes_pred - tes_real), axis=-1)
    val_cor = pearsonr(val_pred, val_real)[0]
    tes_cor = pearsonr(tes_pred, tes_real)[0]
    tra_cor = pearsonr(tra_pred, tra_real)[0]
    tra_mae = np.mean(np.absolute(tra_pred - tra_real), axis=-1) * t_sigma
    val_mae = np.mean(np.absolute(val_pred - val_real), axis=-1) * t_sigma
    tes_mae = np.mean(np.absolute(tes_pred - tes_real), axis=-1) * t_sigma

    return [
        val_los, tes_los, val_cor, tes_cor, tra_cor, tra_mae, val_mae, tes_mae
    ]


def main(args):
    pass


if __name__ == '__main__':
    main()
