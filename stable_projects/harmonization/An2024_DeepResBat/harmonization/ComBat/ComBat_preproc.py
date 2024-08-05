#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import numpy as np
from copy import deepcopy
from config import global_config
from utils.misc import txt2list, csvs2df, save_df
from utils.imputation import forward_filling


def ComBat_preproc(data_path,
                   harm_output_path,
                   covariates,
                   prefix,
                   train_file,
                   val_file,
                   test_files,
                   nb_folds=10):
    """
    Preprocess data for ComBat harmonization
    1. forward filling for missing covariates
    2. merge train_csv and val_csv for fairness
    3. Change DX from [0, 1, 2] to [1, 2, 3]
    4. extract ROI/sites/covariates for ComBat input
    """
    ROIs = txt2list(global_config.ROI_features_path)
    cols = txt2list(global_config.columns_path)
    if 'SIM' in covariates:
        cols += ['SIM']
    train_csv = train_file + '.csv'
    val_csv = val_file + '.csv'
    test_csvs = [test_file + '.csv' for test_file in test_files]
    for fold in range(nb_folds):
        fold_data_path = os.path.join(data_path, str(fold))
        fold_out_path = os.path.join(harm_output_path, str(fold))
        # read original data
        # merge train and val
        train_val_df = csvs2df(fold_data_path, [train_csv, val_csv])
        train_val_df = train_val_df[cols]

        test_dfs = []
        for test_csv in test_csvs:
            temp_test_df = csvs2df(fold_data_path, [test_csv])
            temp_test_df = temp_test_df[cols]
            temp_test_df['DX'] += 1
            test_dfs.append(temp_test_df)
        # forward filling for missing covariates
        if len(covariates) == 0:
            pass
        else:
            train_val_df, test_dfs = ComBat_preproc_ff_covariates(
                train_val_df, test_dfs, covariates)
        subjects = np.unique(train_val_df.RID)
        if '1195' in subjects:
            mask = (train_val_df.RID == '1195') & (
                train_val_df.EXAMDATE == '2011-02-15')
            index = train_val_df[mask].index
            train_val_df.iloc[index, 2] = 0
            train_val_df.iloc[index, 6] = 30
        train_val_df['DX'] += 1
        test_dfs2 = []
        for test_df in test_dfs:
            test_df2 = deepcopy(test_df)
            subjects = np.unique(test_df2.RID)
            if '1195' in subjects:
                mask = (test_df2.RID == '1195') & (
                    test_df2.EXAMDATE == '2011-02-15')
                index = test_df2[mask].index
                test_df2.iloc[index, 2] = 0
                test_df2.iloc[index, 6] = 30
            test_df2['DX'] += 1
            test_dfs2.append(test_df2)
        # extract ROI/sites/covariates
        dfs2save = [train_val_df] + test_dfs2
        save_names = [prefix + '_train_val'] + test_files
        ComBat_preproc_extract_ROI(dfs2save, ROIs, save_names, fold_out_path)
        ComBat_preproc_extract_site(dfs2save, save_names, fold_out_path)
        ComBat_preproc_extract_covariates(dfs2save, covariates, save_names,
                                          fold_out_path)


def ComBat_preproc_ff_covariates(train_df, test_dfs, covariates):
    """
    Run forward filling for missing covariates
    """
    # train
    for cov in covariates:
        train_df = forward_filling(train_df, cov)
    # tests
    filled_test_dfs = []
    for test_df in test_dfs:
        for cov in covariates:
            test_df = forward_filling(test_df, cov)
        filled_test_dfs.append(test_df)

    return train_df, filled_test_dfs


def ComBat_preproc_extract_ROI(dfs, ROIs, save_names, output_path):
    """
    Extract ROI
    """
    assert len(dfs) == len(save_names), 'dfs is not matched with save_names!'
    for df, save_name in zip(dfs, save_names):
        save_path = os.path.join(output_path, save_name + '_ROI.csv')
        save_df(df[ROIs].T, save_path)


def ComBat_preproc_extract_site(dfs, save_names, output_path):
    """
    Extract site
    """
    assert len(dfs) == len(save_names), 'dfs is not matched with save_names!'
    for df, save_name in zip(dfs, save_names):
        save_path = os.path.join(output_path, save_name + '_site.csv')
        save_df(df['SITE'], save_path, ['SITE'])


def ComBat_preproc_extract_covariates(dfs, covariates, save_names,
                                      output_path):
    """
    Extract covariates
    """
    assert len(dfs) == len(save_names), 'dfs is not matched with save_names!'
    for df, save_name in zip(dfs, save_names):
        save_path = os.path.join(output_path, save_name + '_covariates.csv')
        save_df(df[covariates], save_path)
