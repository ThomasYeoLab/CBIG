#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import sys
import random
import argparse
import numpy as np
import pandas as pd
import scipy.io as sio
from sklearn.datasets import make_regression
sys.path.append('../cbig/Chen2024')
from config import config


def save_meta_set(dir_out, eid, phe, fc, df, name):
    name = 'exp_' + name
    os.makedirs(os.path.join(dir_out, name), exist_ok=True)
    np.savetxt(os.path.join(dir_out, name, name + '_subj_list.txt'), eid, fmt='%s')
    np.savetxt(os.path.join(dir_out, name, name + '_phe_list.txt'), phe, fmt='%s')
    sio.savemat(os.path.join(dir_out, name, name + '_fc.mat'), dict([('corr_mat', fc)]))
    df_tmp = df.copy()
    df_tmp = df_tmp[['eid'] + phe]
    df_tmp = df_tmp.loc[df_tmp['eid'].isin(eid)]
    df_tmp.to_csv(os.path.join(dir_out, name, name + '_phe_tab.csv'))


def generate_data(dataset):
    if dataset == 'exp':
        dir_out = config.IN_DIR_EXP
        os.makedirs(dir_out, exist_ok=True)

        # parameter for data generation
        n_samples_XL, n_samples_L, n_samples_M, n_samples_test = 5000, 2000, 1000, 1000
        n_roi = 20
        n_features = int(n_roi * (n_roi - 1) / 2)
        n_informative = 10
        n_targets_XL, n_targets_L, n_targets_M, n_targets_test = 10, 10, 10, 10

    # generate data
    x_XL, y_XL = make_regression(
        n_samples=n_samples_XL,
        n_features=n_features,
        n_informative=n_informative,
        n_targets=n_targets_XL,
        noise=random.random(),
        random_state=seed)

    x_L, y_L = make_regression(
        n_samples=n_samples_L,
        n_features=n_features,
        n_informative=n_informative,
        n_targets=n_targets_L,
        noise=random.random(),
        random_state=seed)

    x_M, y_M = make_regression(
        n_samples=n_samples_M,
        n_features=n_features,
        n_informative=n_informative,
        n_targets=n_targets_M,
        noise=random.random(),
        random_state=seed)

    x_test, y_test = make_regression(
        n_samples=n_samples_test,
        n_features=n_features,
        n_informative=n_informative,
        n_targets=n_targets_test,
        noise=random.random(),
        random_state=seed)

    # convert x to FC format
    index = np.tril(np.ones(n_roi), k=-1) == 1
    fc_list = []
    for i in range(n_samples_XL):
        temp = np.zeros((n_roi, n_roi))
        temp_flat = x_XL[i, :]
        temp[index] = temp_flat
        temp = temp + temp.T
        fc_list.append(temp)
    fc_XL = np.stack(fc_list, axis=2)
    print('converted to fc with shape', fc_XL.shape)
    fc_list = []
    for i in range(n_samples_L):
        temp = np.zeros((n_roi, n_roi))
        temp_flat = x_L[i, :]
        temp[index] = temp_flat
        temp = temp + temp.T
        fc_list.append(temp)
    fc_L = np.stack(fc_list, axis=2)

    fc_list = []
    for i in range(n_samples_M):
        temp = np.zeros((n_roi, n_roi))
        temp_flat = x_M[i, :]
        temp[index] = temp_flat
        temp = temp + temp.T
        fc_list.append(temp)
    fc_M = np.stack(fc_list, axis=2)

    fc_list = []
    for i in range(n_samples_test):
        temp = np.zeros((n_roi, n_roi))
        temp_flat = x_M[i, :]
        temp[index] = temp_flat
        temp = temp + temp.T
        fc_list.append(temp)
    fc_test = np.stack(fc_list, axis=2)

    # generate eid and phenotype name
    eid_XL = np.arange(n_samples_XL) + 100000
    phe_XL = [str(100 + i) + '-888.0' for i in np.arange(n_targets_XL)]
    eid_L = np.arange(n_samples_L) + 200000
    phe_L = [str(200 + i) + '-888.0' for i in np.arange(n_targets_L)]
    eid_M = np.arange(n_samples_M) + 300000
    phe_M = [str(300 + i) + '-888.0' for i in np.arange(n_targets_M)]
    eid_test = np.arange(n_samples_test) + 400000
    phe_test = [str(400 + i) + '-888.0' for i in np.arange(n_targets_test)]

    # put y into dataframe with eid and name
    df_XL = pd.DataFrame(y_XL, columns=phe_XL)
    df_XL.insert(0, "eid", eid_XL)

    df_L = pd.DataFrame(y_L, columns=phe_L)
    df_L.insert(0, "eid", eid_L)

    df_M = pd.DataFrame(y_M, columns=phe_M)
    df_M.insert(0, "eid", eid_M)

    df_test = pd.DataFrame(y_test, columns=phe_test)
    df_test.insert(0, "eid", eid_test)

    # save out data generated
    save_meta_set(dir_out, eid_XL, phe_XL, fc_XL, df_XL, 'train_XL')
    save_meta_set(dir_out, eid_L, phe_L, fc_L, df_L, 'train_L')
    save_meta_set(dir_out, eid_M, phe_M, fc_M, df_M, 'train_M')
    save_meta_set(dir_out, eid_test, phe_test, fc_test, df_test, 'test')

    np.savez(os.path.join(dir_out, 'exp_train_XL', 'exp_train_XL_dnn_input.npz'),
             x_raw=x_XL, y_raw=df_XL[phe_XL].values)
    np.savez(os.path.join(dir_out, 'exp_train_L', 'exp_train_L_dnn_input.npz'),
             x_raw=x_L, y_raw=df_L[phe_L].values)
    np.savez(os.path.join(dir_out, 'exp_train_M', 'exp_train_M_dnn_input.npz'),
             x_raw=x_M, y_raw=df_M[phe_M].values)
    np.savez(os.path.join(dir_out, 'exp_test', 'exp_test_dnn_input.npz'),
             x_raw=x_test, y_raw=df_test[phe_test].values)


if __name__ == "__main__":
    seed = config.RAMDOM_SEED
    random.seed(seed)
    np.random.seed(seed)

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset", type=str, default="exp")
    args = parser.parse_args()

    generate_data(args.dataset)
