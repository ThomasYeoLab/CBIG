#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import warnings
import numpy as np
import pandas as pd
from config import global_config
from utils.misc import save_df, create_folder

warnings.filterwarnings('ignore')


def extract_data_from_split(data_df, split_df, isMatched, isTrain):
    """
    Extract data using split_df

    Args:
        data_df (class DataFrame): Data
        split_df (class DataFrame): Split information
        isMatched (bool): Whether the subject is matched
        isTrain (bool): Whether the subject is for training
    """
    split_df_mask = (split_df.isMatched == isMatched) & (
        split_df.isTrain == isTrain)
    split_df = split_df[split_df_mask]
    split_df['RID'] = split_df['RID'].apply(str)
    data_df['RID'] = data_df['RID'].apply(str)
    cols = list(data_df.columns)
    val_df = pd.DataFrame(columns=cols)
    for i in range(split_df.shape[0]):
        rid = split_df.iloc[i, 0]
        exam_date = split_df.iloc[i, 1]
        mask = (data_df.RID == rid) & (data_df.EXAMDATE == exam_date)
        val_df = val_df.append(data_df[mask], ignore_index=True)
    return val_df


def data_split(data_path, dataset_pair, nb_folds=10):
    """
    Split data into <nb_folds> using given split config file

    Args:
        data_path (str): Path for data
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
        nb_folds (int, optional): Number of folds. Defaults to 10.
    """
    # output_path
    output_path = os.path.join(data_path, 'splits', dataset_pair)
    # read matched and unmatched data
    nonADNI_name = dataset_pair.split('-')[1]
    ADNI_matched_df = pd.read_csv(
        os.path.join(output_path, 'ADNI_matched.csv'))
    ADNI_unmatched_df = pd.read_csv(
        os.path.join(output_path, 'ADNI_unmatched.csv'))
    nonADNI_matched_df = pd.read_csv(
        os.path.join(output_path, nonADNI_name + '_matched.csv'))
    nonADNI_unmatched_df = pd.read_csv(
        os.path.join(output_path, nonADNI_name + '_unmatched.csv'))
    matched_df = ADNI_matched_df.append(nonADNI_matched_df, ignore_index=True)
    matched_df['RID'] = matched_df['RID'].apply(str)
    unmatched_df = ADNI_unmatched_df.append(
        nonADNI_unmatched_df, ignore_index=True)
    unmatched_df['RID'] = unmatched_df['RID'].apply(str)
    # split val/test for 10 folds
    for fold in range(nb_folds):
        fold_output_path = os.path.join(output_path, str(fold))
        create_folder(fold_output_path)
        # load split config
        split_df = pd.read_csv(
            os.path.join(fold_output_path,
                         dataset_pair + '_fold_' + str(fold) + '_splits.csv'))
        # save testset
        # for unmatch2match, the testset is matched_df
        save_df(matched_df,
                os.path.join(fold_output_path, 'unmatch2match_test.csv'))
        # for match2unmatch, the testset is unmatched_df
        save_df(unmatched_df,
                os.path.join(fold_output_path, 'match2unmatch_test.csv'))
        # extract validation set
        matched_val_df = extract_data_from_split(matched_df, split_df, 1, 0)
        unmatched_val_df = extract_data_from_split(unmatched_df, split_df, 0,
                                                   0)
        # edge case for subject 1195
        if '1195' in list(unmatched_val_df.RID):
            mask = (unmatched_val_df.RID == '1195') & (
                unmatched_val_df.EXAMDATE == '2011-02-15') & (np.isnan(
                    unmatched_val_df.DX))
            index = list(unmatched_val_df[mask].index)[1:]
            unmatched_val_df.drop(index, inplace=True)
            mask = (unmatched_val_df.RID == '1195') & (
                unmatched_val_df.EXAMDATE == '2011-02-15') & (
                    ~np.isnan(unmatched_val_df.DX))
            index = list(unmatched_val_df[mask].index)[1:]
            unmatched_val_df.drop(index, inplace=True)
        # save
        save_df(matched_val_df,
                os.path.join(fold_output_path, 'match2unmatch_val.csv'))
        save_df(unmatched_val_df,
                os.path.join(fold_output_path, 'unmatch2match_val.csv'))
        # extract training data
        matched_train_df = extract_data_from_split(matched_df, split_df, 1, 1)
        unmatched_train_df = extract_data_from_split(unmatched_df, split_df, 0,
                                                     1)
        # edge case for subject 1195
        if '1195' in list(unmatched_train_df.RID):
            mask = (unmatched_train_df.RID == '1195') & (
                unmatched_train_df.EXAMDATE == '2011-02-15') & (np.isnan(
                    unmatched_train_df.DX))
            index = list(unmatched_train_df[mask].index)[1:]
            unmatched_train_df.drop(index, inplace=True)
            mask = (unmatched_train_df.RID == '1195') & (
                unmatched_train_df.EXAMDATE == '2011-02-15') & (
                    ~np.isnan(unmatched_train_df.DX))
            index = list(unmatched_train_df[mask].index)[:1]
            unmatched_train_df.drop(index, inplace=True)
        save_df(unmatched_train_df,
                os.path.join(fold_output_path, 'unmatch2match_train.csv'))
        save_df(matched_train_df,
                os.path.join(fold_output_path, 'match2unmatch_train.csv'))


if __name__ == '__main__':
    data_path = os.path.join(global_config.root_path, 'data')
    data_split(data_path, 'ADNI-AIBL')
    data_split(data_path, 'ADNI-MACC')
