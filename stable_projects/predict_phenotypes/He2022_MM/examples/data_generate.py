#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import random
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


def generate_data(dir_out):
    # parameter for data generation
    n_samples = 1000
    n_roi = 10
    n_features = int(n_roi * (n_roi - 1) / 2)
    n_informative = int(n_features / 2)
    n_targets = 10

    # generate data
    x, y = make_regression(
        n_samples=n_samples,
        n_features=n_features,
        n_informative=n_informative,
        n_targets=n_targets,
        noise=0.1,
        random_state=seed)

    print('randomly generated x with shape', x.shape, 'y with shape', y.shape)
    print('x (part of it):\n', x[:10, :5])
    print('y (part of it):\n', y[:10, :5])

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
    ind_tra = range(0, int(n_samples / 2))
    ind_tes = range(int(n_samples / 2), n_samples)
    inds = range(0, n_samples)
    eid_tra = eid[ind_tra]
    eid_tes = eid[ind_tes]
    phe_tra = phes[:int(n_targets / 2)]
    phe_tes = phes[int(n_targets / 2):]

    # save out data generated
    save_meta_set(dir_out, eid_tes, phe_tes, ind_tes, pfc, df, 'test')
    save_meta_set(dir_out, eid_tra, phe_tra, ind_tra, pfc, df, 'train')
    save_meta_set(dir_out, eid, phes, inds, pfc, df, 'train_test')


if __name__ == "__main__":
    seed = 0
    random.seed(seed)
    np.random.seed(seed)

    base_dir = os.path.join(
        os.getenv('CBIG_CODE_DIR'),
        'stable_projects/predict_phenotypes/He2022_MM')
    dir_out = os.path.join(base_dir, 'examples/exp_input')
    os.makedirs(dir_out, exist_ok=True)
    generate_data(dir_out)
