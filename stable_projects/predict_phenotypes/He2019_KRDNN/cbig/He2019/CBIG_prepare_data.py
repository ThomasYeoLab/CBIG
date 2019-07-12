#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import numpy as np
import scipy.io as sio
from .config import config


def data_ukbb_base(corr_mat=None,
                   beh=None,
                   is_sex=False,
                   should_flatten=False,
                   is_gcnn=False):
    '''data load and process base function for UK Biobank dataset

    Args:
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.
        is_sex (bool, optional): whether data y is sex
        should_flatten (bool, optional): whether flatten data X (Functional
            connectivity)
        is_gcnn (bool, optional): whether return value for GCNN

    Returns:
        Tuple: if is_gcnn set to False, tuple returned contains splitted x
            and y data for training, validation and testing. If is_gcnn set
            to True, tuple returned contains x and y data and index (mask) for
            training, validation and testing.
    '''
    if corr_mat is None:
        n_subjects = config.UKBB_NUM_SUBJECT
        measure_sets = config.UKBB_MEASURE_SETS
        mat_content = sio.loadmat(
            os.path.join(config.UKBB_ORIG_DIR, config.UKBB_CORR_MAT))
        x = np.transpose(mat_content['corr_mat'], (2, 0, 1))

        # load sub-fold info
        mat_content = sio.loadmat(
            os.path.join(config.UKBB_ORIG_DIR, config.UKBB_SUBJECT_LIST))
        subject_split = mat_content['subject_split']
        index_tra = np.squeeze(subject_split['index_tra'][0][0]).astype(bool)
        index_val = np.squeeze(subject_split['index_val'][0][0]).astype(bool)
        index_tes = np.squeeze(subject_split['index_tes'][0][0]).astype(bool)
    else:
        n_subjects = config.EXAMPLE_N_SUBJECT
        measure_sets = ['nilearn_adhd']
        x = corr_mat

        # sub-fold
        p1 = int(n_subjects / 4)
        p2 = int(n_subjects / 2)
        p3 = p1 + p2
        index_tra = np.concatenate((np.ones((p2)), np.zeros(
            (p2)))).astype(bool)
        index_val = np.concatenate((np.zeros((p2)), np.ones((p1)),
                                    np.zeros((p1)))).astype(bool)
        index_tes = np.concatenate((np.zeros((p3)), np.ones(
            (p1)))).astype(bool)

    # flatten and index x
    if should_flatten:
        tmp = np.tril_indices(x.shape[1], -1)
        x = x[:, tmp[0], tmp[1]]
        train_x = x[index_tra, :]
        valid_x = x[index_val, :]
        test_x = x[index_tes, :]
    else:
        train_x = np.expand_dims(x[index_tra, :, :], axis=1)
        valid_x = np.expand_dims(x[index_val, :, :], axis=1)
        test_x = np.expand_dims(x[index_tes, :, :], axis=1)

    # load y
    temp_concat = np.zeros((n_subjects, 0))
    for j in range(len(measure_sets)):
        if beh is None:
            mat_content = sio.loadmat(
                os.path.join(config.UKBB_ORIG_DIR, 'y', 'dis_y_regress',
                             measure_sets[j] + '_y_regress.mat'))
            beh = mat_content['measures_num_regress']
        if is_sex:
            temp = np.expand_dims(beh[:, 0], -1)
            assert np.unique(temp).shape[0] == 2, "More than 2 class."
            temp = temp - np.min(temp)
        else:
            temp = beh
        temp_concat = np.concatenate((temp_concat, temp), axis=1)
    y = temp_concat

    # index y
    train_y = y[index_tra, :]
    valid_y = y[index_val, :]
    test_y = y[index_tes, :]

    if is_gcnn:
        return x, y, index_tra, index_val, index_tes
    else:
        return train_x, train_y, valid_x, valid_y, test_x, test_y


def data_ukbb_fnn(save_dir, corr_mat=None, beh=None):
    '''data load and process function for FNN and UK Biobank dataset

    Args:
        save_dir (str): directory to save the data
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.

    Returns:
        None
    '''
    train_x, train_y, valid_x, valid_y, test_x, test_y = data_ukbb_base(
        corr_mat, beh, should_flatten=True)
    np.savez(
        os.path.join(save_dir, 'data_fnn.npz'),
        test_x=test_x,
        test_y=test_y,
        train_x=train_x,
        train_y=train_y,
        valid_x=valid_x,
        valid_y=valid_y)


def data_ukbb_fnn_sex(save_dir, corr_mat=None, beh=None):
    '''data load and process function for FNN, sex and UK Biobank dataset

    Args:
        save_dir (str): directory to save the data
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.

    Returns:
        None
    '''
    train_x, train_y, valid_x, valid_y, test_x, test_y = data_ukbb_base(
        corr_mat, beh, is_sex=True, should_flatten=True)
    np.savez(
        os.path.join(save_dir, 'data_fnn_sex.npz'),
        test_x=test_x,
        test_y=test_y,
        train_x=train_x,
        train_y=train_y,
        valid_x=valid_x,
        valid_y=valid_y)


def data_ukbb_brainnetcnn(save_dir, corr_mat=None, beh=None):
    '''data load and process function for BrainNetCNN and UK Biobank dataset

    Args:
        save_dir (str): directory to save the data
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.

    Returns:
        None
    '''
    train_x, train_y, valid_x, valid_y, test_x, test_y = data_ukbb_base(
        corr_mat, beh)
    np.savez(
        os.path.join(save_dir, 'data_brainnetcnn.npz'),
        test_x=test_x,
        test_y=test_y,
        train_x=train_x,
        train_y=train_y,
        valid_x=valid_x,
        valid_y=valid_y)


def data_ukbb_brainnetcnn_sex(save_dir, corr_mat=None, beh=None):
    '''data load and process function for BrainNetCNN, sex and UK Biobank dataset

    Args:
        save_dir (str): directory to save the data
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.

    Returns:
        None
    '''
    train_x, train_y, valid_x, valid_y, test_x, test_y = data_ukbb_base(
        corr_mat, beh, is_sex=True)
    np.savez(
        os.path.join(save_dir, 'data_brainnetcnn_sex.npz'),
        test_x=test_x,
        test_y=test_y,
        train_x=train_x,
        train_y=train_y,
        valid_x=valid_x,
        valid_y=valid_y)


def data_ukbb_gcnn(save_dir, corr_mat=None, beh=None):
    '''data load and process function for GCNN and UK Biobank dataset

    Args:
        save_dir (str): directory to save the data
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.

    Returns:
        None
    '''
    x, y, index_tra, index_val, index_tes = data_ukbb_base(
        corr_mat, beh, should_flatten=True, is_gcnn=True)
    np.savez(
        os.path.join(save_dir, 'data_gcnn.npz'),
        input_x=x,
        input_y=y,
        test_mask=index_tes,
        train_mask=index_tra,
        valid_mask=index_val)


def data_ukbb_gcnn_sex(save_dir, corr_mat=None, beh=None):
    '''data load and process function for GCNN and UK Biobank dataset

    Args:
        save_dir (str): directory to save the data
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.

    Returns:
        None
    '''
    x, y, index_tra, index_val, index_tes = data_ukbb_base(
        corr_mat, beh, is_sex=True, should_flatten=True, is_gcnn=True)
    np.savez(
        os.path.join(save_dir, 'data_gcnn_sex.npz'),
        input_x=x,
        input_y=y,
        test_mask=index_tes,
        train_mask=index_tra,
        valid_mask=index_val)


def data_hcp_base(n_folds, corr_mat=None, beh=None, should_flatten=False):
    '''data load and process base function for HCP dataset

    Args:
        n_folds (int): number of folds
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.
        should_flatten (bool, optional): whether flatten data X (Functional
            connectivity)

    Returns:
        Tuple: contains data x, list of index and y data for cross validation
            folds
    '''
    if corr_mat is None:
        n_subjects = config.HCP_NUM_SUBJECT
        # load x (correlation matrix)
        mat_content = sio.loadmat(
            os.path.join(config.HCP_ORIG_DIR, config.HCP_CORR_MAT))
        x = np.transpose(mat_content['corr_mat'], (2, 0, 1))

        # load sub-fold info
        mat_content = sio.loadmat(
            os.path.join(config.HCP_ORIG_DIR, config.HCP_SUBJECT_LIST))
        sub_fold = mat_content['sub_fold']
    else:
        n_subjects = config.EXAMPLE_N_SUBJECT
        x = corr_mat

    list_fold_index = []
    list_fold_y = []
    for i in range(n_folds):
        if corr_mat is None:
            # load sub-fold list
            tmp = np.concatenate(sub_fold[i]['fold_index'][0]).astype(int)
            list_fold_index.append(tmp)
            # load y for this sub-fold
            tmp = os.path.join(config.HCP_ORIG_DIR, 'y', 'fold_' + str(i + 1),
                               'y_regress.mat')
            mat_content = sio.loadmat(tmp)
            fold_y = mat_content['y_resid']
            list_fold_y.append(fold_y)
        else:
            temp = np.zeros((n_subjects))
            cnt = int(n_subjects / n_folds)
            temp[(i * cnt):(i * cnt + cnt)] = 1
            list_fold_index.append(temp.astype(bool))
            list_fold_y.append(beh)

    # flatten x
    if should_flatten:
        tmp = np.tril_indices(x.shape[1], -1)
        x = x[:, tmp[0], tmp[1]]

    return x, list_fold_index, list_fold_y


def data_hcp_fnn(save_dir, n_folds, corr_mat=None, beh=None):
    '''data load and save function for HCP dataset and FNN network

    Args:
        save_dir (str): directory to save the data
        n_folds (int): number of folds
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.

    Returns:
        None
    '''
    x, index, y = data_hcp_base(n_folds, corr_mat, beh, should_flatten=True)
    os.makedirs(save_dir + '/fnn', exist_ok=True)
    for i in range(n_folds):
        # index test data for this fold
        test_index = np.nonzero(index[i])[0]
        test_x = x[test_index, :]
        test_y = y[i][test_index, :]
        # index train and validation fold
        train_valid_x = []
        train_valid_y = []
        for j in range(n_folds):
            if j == i:
                continue
            train_valid_index = np.nonzero(index[j])[0]
            train_valid_x.append(x[train_valid_index, :])
            train_valid_y.append(y[i][train_valid_index, :])
        # save data
        np.savez(
            os.path.join(save_dir, 'fnn', 'data_fold' + str(i + 1) + '.npz'),
            test_x=test_x,
            test_y=test_y,
            train_valid_x=train_valid_x,
            train_valid_y=train_valid_y)
    return


def data_hcp_brainnetcnn(save_dir, n_folds, corr_mat=None, beh=None):
    '''data load and save function for HCP dataset and BrainNetCNN network

    Args:
        save_dir (str): directory to save the data
        n_folds (int): number of folds
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.

    Returns:
        None
    '''
    x, index, y = data_hcp_base(n_folds, corr_mat, beh)
    os.makedirs(save_dir + '/brainnetcnn', exist_ok=True)
    for i in range(n_folds):
        # index test data for this fold
        test_index = np.nonzero(index[i])[0]
        test_x = x[test_index, :, :]
        test_y = y[i][test_index, :]
        # index train and validation fold
        train_valid_x = []
        train_valid_y = []
        for j in range(n_folds):
            if j == i:
                continue
            train_valid_index = np.nonzero(index[j])[0]
            train_valid_x.append(x[train_valid_index, :, :])
            train_valid_y.append(y[i][train_valid_index, :])
        # save data
        np.savez(
            os.path.join(save_dir, 'brainnetcnn',
                         'data_fold' + str(i + 1) + '.npz'),
            test_x=test_x,
            test_y=test_y,
            train_valid_x=train_valid_x,
            train_valid_y=train_valid_y)
    return


def data_hcp_gcnn(save_dir, n_folds, corr_mat=None, beh=None):
    '''data load and save function for HCP dataset and GCNN network

    Args:
        save_dir (str): directory to save the data
        n_folds (int): number of folds
        corr_mat (ndarray, optional): functional connectivity matrix for all
            subjects. Only used for example data.
        beh (ndarray, optional): behavioral measures for all
            subjects. Only used for example data.

    Returns:
        None
    '''
    x, index, y = data_hcp_base(n_folds, corr_mat, beh, should_flatten=True)
    os.makedirs(save_dir + '/gcnn', exist_ok=True)
    for i in range(n_folds):
        temp = index.copy()
        test_mask = temp.pop(i)
        train_valid_mask = temp
        input_x = x
        input_y = y[i]
        # save data
        np.savez(
            os.path.join(save_dir, 'gcnn', 'data_fold' + str(i + 1) + '.npz'),
            input_x=input_x,
            input_y=input_y,
            test_mask=test_mask,
            train_valid_mask=train_valid_mask)

    return


def get_gcnn_graph(graph_dir, x, k=3, d=5):
    """generate adjacency matrix for GCNN

    Args:
        graph_dir (str): directory to save the matrix
        x (ndarray): functional connectivity matrix for all subjects
        k (int, optional): options for generation
            k=1: for each nodes, only keep edges with top d(int) correlation
            k=2: sort the correlation, only keep top d(int)% edges
            k=3: combination of option 1 and option 2, using same d(int) value
        d (int, optional): threshold for generation

    Returns:
        None
    """
    # compute adjancecy matrix
    tmp = np.tril_indices(x.shape[1], -1)
    x = x[:, tmp[0], tmp[1]]
    fsm = np.corrcoef(x)
    n_node = fsm.shape[0]
    fsm[np.diag_indices(n_node)] = 0

    # process adj matrix
    if k == 1:
        for i in range(n_node):
            index = np.argsort(fsm[i, :])[:-d]
            fsm[i, index] = 0
    elif k == 2:
        index = int(np.floor(d * n_node * n_node / 100) + 1)
        index = np.unravel_index(
            np.argsort(fsm, axis=None)[:-index], fsm.shape)
        fsm[index] = 0
    elif k == 3:
        mark = np.zeros((n_node, n_node))
        index = int(np.floor(d * n_node * n_node / 100) + 1)
        index = np.unravel_index(
            np.argsort(fsm, axis=None)[-index:], fsm.shape)
        mark[index] = 1
        for i in range(n_node):
            index = np.argsort(fsm[i, :])[-d:]
            mark[i, index] = 1
        fsm = np.multiply(fsm, mark)
    else:
        assert False, "wrong graph option"

    # undirected graph
    tmp = fsm.T > fsm
    fsm = fsm - np.multiply(fsm, tmp) + np.multiply(fsm.T, tmp)

    # analyze the graph
    cnt = np.sum(fsm > 0, axis=0)
    print('vertex degree mean:', np.mean(cnt), 'min:', np.min(cnt), 'max:',
          np.max(cnt))

    # check connectivity
    adj = fsm
    while 1:
        adj_o = adj
        adj = np.logical_or(adj, np.matmul(adj, adj.T)).astype(float)
        if np.min((adj_o == adj).astype(float)) == 1:
            break
    if np.min(adj) == 0:
        assert False, 'Error: not fully connected'
    else:
        print('fully connected')

    # save graph
    os.makedirs(graph_dir, exist_ok=True)
    file = os.path.join(
        graph_dir,
        str(n_node) + 'Subject_corr_option_' + str(k) + '_param_' + str(d) +
        '.cites')
    f = open(file, 'w')
    rows, cols = np.nonzero(fsm)
    for i in range(rows.shape[0]):
        f.write(
            str(rows[i] + 1) + '\t' + str(cols[i] + 1) + '\t' +
            str(fsm[rows[i], cols[i]]) + '\n')
    f.close()
    return


def main():
    pass


if __name__ == '__main__':
    main()
