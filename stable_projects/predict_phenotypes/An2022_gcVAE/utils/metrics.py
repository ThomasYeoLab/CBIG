#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
from sklearn import metrics


def mae(pred, target, device='cuda'):
    """
    Calculate MAE for 2 tensors (GPU or CPU)

    Args:
        pred (array): Tensor array
        target (array): Tensor array
        device (str, optional): CPU or GPU. Defaults to 'cuda'.
    """
    if device == 'cuda':
        # GPU tensors
        pred = pred.data.cpu().numpy()
        target = target.data.cpu().numpy()
    else:
        # CPU tensors
        pred = pred.numpy()
        target = target.numpy()
    return metrics.mean_absolute_error(pred, target)


def subject_mae(pred, target, RID_column, device='cuda'):
    """
    Compute subject level MAE for 2 tensors

    Args:
        pred (array): Tensor array
        target (array): Tensor array
        RID_column (class Series): Series for RIDs
        device (str, optional): CPU or GPU. Defaults to 'cuda'.
    """
    if device == 'cuda':
        # GPU tensors
        pred = pred.data.cpu().numpy()
        target = target.data.cpu().numpy()
    elif device == 'cpu':
        # CPU tensors
        pred = pred.numpy()
        target = target.numpy()
    else:
        pass
    # calculate subject level MAE
    subs, ind = np.unique(RID_column, return_index=True)
    subs = subs[np.argsort(ind)]
    nb_subs = subs.shape[0]
    # Copy RID_column along axis 1 to get a #TPsX#Subs matrix
    RID_matrix = np.tile(RID_column, (1, nb_subs))
    mask = (~(RID_matrix == subs)) * 1  # 0 means select, 1 means not select
    # copy pred along axis 1, generate a #TPsX#Subs matrix
    pred_matrix = np.tile(pred, (1, nb_subs))
    # get absolute error matrix
    ae_matrix = abs(pred_matrix - target)
    mae_vector = np.mean(np.ma.masked_array(ae_matrix, mask), axis=0)
    return mae_vector, np.mean(mae_vector)


def subject_acc(pred_prob, target, RID_column, device='cuda'):
    """
    Compute subject level accuracy for 2 tensors

    Args:
        pred_prob (array): Tensor array
        target (array): Tensor array
        RID_column (class Series): Series for RIDs
        device (str, optional): CPU or GPU. Defaults to 'cuda'.
    """
    if device == 'cuda':
        # GPU tensors
        pred_prob = pred_prob.data.cpu().numpy()
        target = target.data.cpu().numpy()
    elif device == 'cpu':
        # CPU tensors
        pred_prob = pred_prob.numpy()
        target = target.numpy()
    else:
        # numpy arrays
        pass
    # calculate subject level accuracy
    pred = pred_prob.argmax(axis=1)
    pred = pred.reshape((pred.shape[0], 1))
    subs, ind = np.unique(RID_column, return_index=True)
    subs = subs[np.argsort(ind)]
    nb_subs = subs.shape[0]
    # Copy RID_column along axis 1 to get a #TPsX#Subs matrix
    RID_matrix = np.tile(RID_column, (1, nb_subs))
    mask = (~(RID_matrix == subs)) * 1  # 0 means select, 1 means not select
    # copy pred along axis 1, generate a #TPsX#Subs matrix
    pred_matrix = np.tile(pred, (1, nb_subs))
    # get absolute error matrix
    acc_matrix = ((pred_matrix == target) * 1.0)
    acc_vector = np.mean(np.ma.masked_array(acc_matrix, mask), axis=0)

    return acc_vector, np.mean(acc_vector)


def site_prediction_metric(RID_column, pred_prob, truth, threshold):
    """
    Compute subject level site prediction accuracy

    Args:
        RID_column (class Series): Series for RIDs
        pred_prob (ndarray): Numpy array
        truth (ndarray): Numpy array
        threshold (ndarray): Numpy array
    """
    # get prediction
    pred = (pred_prob >= threshold) * 1.
    subs = np.unique(RID_column)
    subs, ind = np.unique(RID_column, return_index=True)
    subs = subs[np.argsort(ind)]
    nb_subs = subs.shape[0]
    # Copy RID_column along axis 1, generate a matrix #Tps by #Subs
    RID_matrix = np.tile(RID_column, (nb_subs, 1)).transpose()
    mask = (~(RID_matrix == subs)) * 1.  # 0 means select, 1 means not select
    # copy pred along axis 1, generate a matrix #Tps by #Subs
    pred_matrix = np.tile(pred, (nb_subs, 1)).transpose()
    # 0 means wrong prediction, 1 means correct prediction
    pred_acc_matrix = ((pred_matrix.transpose() == truth) * 1.0).transpose()
    acc_vector = np.mean(np.ma.masked_array(pred_acc_matrix, mask), axis=0)
    acc_vector = np.array(acc_vector)
    return acc_vector, np.mean(acc_vector)


def site_prediction_acc_auc(pred_prob, truth, threshold):
    """
    Compute timepoint level site prediction accuracy & AUC

    Args:
        pred_prob (ndarray): Numpy array
        truth (ndarray): Numpy array
        threshold (float): Threshold
    """
    pred = (pred_prob >= threshold) * 1
    acc = metrics.accuracy_score(truth, pred)
    auc = metrics.roc_auc_score(truth, pred_prob)

    return acc, auc
