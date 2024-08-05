#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
from config import global_config
from utils.misc import txt2list, csvs2df, \
    replace_with_harmed_ROI, save_df


def ComBat_postproc(data_path,
                    output_path,
                    harm_files,
                    nb_folds=10,
                    isSIM=False):
    """
    Postprocess harmonized data
    """
    ROIs = txt2list(global_config.ROI_features_path)
    features = ['RID', 'EXAMDATE', 'DX', 'SITE', 'AGE', 'SEX', 'MMSE'] + ROIs
    if isSIM:
        features = ['RID', 'EXAMDATE', 'DX', 'SITE', 'AGE', 'SEX', 'SIM'
                    ] + ROIs
    for fold in range(nb_folds):
        fold_data_path = os.path.join(data_path, str(fold))
        fold_output_path = os.path.join(output_path, str(fold))
        for harm_file in harm_files:
            harm_csvs = [harm_file + '.csv']
            origin_csvs = []
            for origin_file in global_config.ComBat_harm_mapper[harm_file]:
                origin_csvs.append(origin_file + '.csv')
            ComBat_potproc_single_csv(fold_data_path, fold_output_path,
                                      features, harm_csvs, origin_csvs)


def ComBat_potproc_single_csv(fold_data_path, fold_output_path, features,
                              harm_csvs, origin_csvs):
    """
    Post process ComBat harmonized data
    """
    if len(origin_csvs) == 2:
        # train_val
        harm_train_val_ROI_df = csvs2df(fold_output_path,
                                        harm_csvs,
                                        index_col=0)
        train_df = csvs2df(fold_data_path, origin_csvs[:1])
        val_df = csvs2df(fold_data_path, origin_csvs[1:])
        # replace with harmonized ROI
        # for train
        harm_train_df = replace_with_harmed_ROI(
            train_df,
            harm_train_val_ROI_df.iloc[:, :train_df.shape[0]],
            features,
            isTranspose=True)
        save_df(harm_train_df, os.path.join(fold_output_path, origin_csvs[0]))
        # for val
        harm_val_df = replace_with_harmed_ROI(
            val_df,
            harm_train_val_ROI_df.iloc[:, train_df.shape[0]:],
            features,
            isTranspose=True)
        save_df(harm_val_df, os.path.join(fold_output_path, origin_csvs[1]))
    else:
        # test
        harm_test_ROI_df = csvs2df(fold_output_path, harm_csvs, index_col=0)
        test_df = csvs2df(fold_data_path, origin_csvs)
        harm_test_df = replace_with_harmed_ROI(test_df,
                                               harm_test_ROI_df,
                                               features,
                                               isTranspose=True)
        save_df(harm_test_df, os.path.join(fold_output_path, origin_csvs[0]))
