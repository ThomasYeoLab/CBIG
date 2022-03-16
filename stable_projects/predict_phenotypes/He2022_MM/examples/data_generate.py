#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import random
import argparse
import numpy as np
import pandas as pd
import scipy.io as sio
from sklearn.datasets import make_regression


def save_meta_set(dir_out, eid, phe, ind, pfc, df, name):
    name = 'exp_' + name
    np.savetxt(os.path.join(dir_out, name + '_subj_list.txt'), eid, fmt='%s')
    np.savetxt(
        os.path.join(dir_out, name + '_final_phe_list.txt'), phe, fmt='%s')
    sio.savemat(
        os.path.join(dir_out, name + '_pfc.mat'),
        dict([('corr_mat', pfc[:, :, ind])]))
    df_tmp = df.copy()
    df_tmp = df_tmp[['eid'] + phe]
    df_tmp = df_tmp.loc[df_tmp['eid'].isin(eid)]
    df_tmp.to_csv(os.path.join(dir_out, name + '_final.csv'))
    return df_tmp


def generate_data(dataset):
    base_dir = os.path.join(
        os.getenv('CBIG_CODE_DIR'),
        'stable_projects/predict_phenotypes/He2022_MM')
    if dataset == 'exp':
        dir_out = os.path.join(base_dir, 'examples/exp_input')
        os.makedirs(dir_out, exist_ok=True)

        # parameter for data generation
        n_samples = 5000
        n_tra_samples = 4000
        n_roi = 20
        n_features = int(n_roi * (n_roi - 1) / 2)
        n_informative = 40
        n_targets = 50
        n_tra_targets = 40
    elif dataset == 'unit_tests':
        dir_out = os.path.join(base_dir, 'examples/unit_tests_input')
        os.makedirs(dir_out, exist_ok=True)

        # parameter for data generation
        n_samples = 800
        n_tra_samples = 400
        n_roi = 8
        n_features = int(n_roi * (n_roi - 1) / 2)
        n_informative = 8
        n_targets = 5
        n_tra_targets = 3

    # generate data
    x, y = make_regression(
        n_samples=n_samples,
        n_features=n_features,
        n_informative=n_informative,
        n_targets=n_targets,
        noise=0.1,
        random_state=seed)

    print('randomly generated x with shape', x.shape, 'y with shape', y.shape)

    # convert x to FC format
    index = np.tril(np.ones(n_roi), k=-1) == 1
    pfc_list = []
    for i in range(n_samples):
        temp = np.zeros((n_roi, n_roi))
        temp_flat = x[i, :]
        temp[index] = temp_flat
        temp = temp + temp.T
        pfc_list.append(temp)
    pfc = np.stack(pfc_list, axis=2)
    print('converted to fc with shape', pfc.shape)

    # generate eid and phenotype name
    eid = np.arange(n_samples) + 8000001
    phes = [str(101 + i) + '-888.0' for i in np.arange(n_targets)]
    # print('generated eid and phe \n', eid, '\n', phes)

    # put y into dataframe with eid and name
    df = pd.DataFrame(y, columns=phes)
    df['eid'] = eid
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    print(df)

    # split data and phenotype
    ind_tra = range(0, n_tra_samples)
    ind_tes = range(n_tra_samples, n_samples)
    inds = range(0, n_samples)
    eid_tra = eid[ind_tra]
    eid_tes = eid[ind_tes]
    phe_tra = phes[:n_tra_targets]
    phe_tes = phes[n_tra_targets:]

    # save out data generated
    df_tes = save_meta_set(dir_out, eid_tes, phe_tes, ind_tes, pfc, df, 'test')
    _ = save_meta_set(dir_out, eid_tra, phe_tra, ind_tra, pfc, df, 'train')
    save_meta_set(dir_out, eid, phes, inds, pfc, df, 'train_test')

    np.savez(
        os.path.join(dir_out, 'exp_dnn_input_test.npz'),
        x_train_raw=x,
        y_train_raw=df[phe_tra].values,
        x_test=x[ind_tes, :],
        y_test_different_set_phe=df_tes[phe_tes].values,
        tra_phe=phe_tra,
        tes_phe=phe_tes)


if __name__ == "__main__":
    seed = 0
    random.seed(seed)
    np.random.seed(seed)

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset", type=str, default="exp")
    args = parser.parse_args()

    generate_data(args.dataset)
