#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import numpy as np
import pandas as pd
from copy import deepcopy
from utils.misc import create_folder, txt2list


def merge_multi_bins(args,
                     round1_threshold_list,
                     round2_threshold_list,
                     name1='MACC',
                     name2='ADNI'):
    """
    Merge multiple bins data

    Args:
        args (tuple): Parameters
        round1_threshold_list (list): List of thresholds in round 1 matching
        round2_threshold_list (list): List of thresholds in round 1 matching
        name1 (str, optional): Name for dataset1. Defaults to 'MACC'.
        name2 (str, optional): Name for dataset2. Defaults to 'ADNI'.
    """
    bins_path = os.path.join(args.checkpoint_path, args.matching_pair,
                             'matching_' + str(args.nb_bins) + 'BINs')
    merged_path = os.path.join(
        bins_path, 'MERGED',
        'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
        + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty))
    create_folder(merged_path)
    controled_MACC = pd.DataFrame()
    controled_ADNI = pd.DataFrame()
    for bin in range(args.nb_bins):
        threhsold_name = 'Threshold_' + str(
            round1_threshold_list[bin]) + '-' + str(round2_threshold_list[bin])
        bin_path = os.path.join(bins_path, 'BIN_' + str(bin), 'merged')
        bin_merged_path = os.path.join(
            bin_path,
            'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) +
            '_DX' + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty),
            threhsold_name)
        bin_macc_df = pd.read_csv(
            os.path.join(bin_merged_path, 'controled_' + name1 + '.csv'))
        bin_adni_df = pd.read_csv(
            os.path.join(bin_merged_path, 'controled_' + name2 + '.csv'))
        controled_MACC = controled_MACC.append(bin_macc_df, sort=False)
        controled_ADNI = controled_ADNI.append(bin_adni_df, sort=False)
    # save merged
    controled_MACC.to_csv(
        os.path.join(merged_path, 'controled_' + name1 + '.csv'),
        sep=',',
        index=False)
    controled_ADNI.to_csv(
        os.path.join(merged_path, 'controled_' + name2 + '.csv'),
        sep=',',
        index=False)


def extract_matched_unmatched_data(args,
                                   input_path,
                                   name1='MACC',
                                   name2='ADNI'):
    """
    Extract matched and unmatched data

    Args:
        args (tuple): Parameters
        input_path (str): Path for input data
        name1 (str, optional): Name for dataset1. Defaults to 'MACC'.
        name2 (str, optional): Name for dataset2. Defaults to 'ADNI'.
    """
    # columns
    cols = txt2list(args.columns_path)
    # read raw data
    raw_dataset1 = pd.read_csv(args.MACC_data_path, usecols=cols)
    raw_dataset2 = pd.read_csv(args.ADNI_data_path, usecols=cols)
    # read matched data
    matched_dataset1 = pd.read_csv(
        os.path.join(input_path, 'controled_' + name1 + '.csv'))
    matched_dataset2 = pd.read_csv(
        os.path.join(input_path, 'controled_' + name2 + '.csv'))
    # deepcopy data
    matched_dataset1_copy = deepcopy(matched_dataset1)
    matched_dataset2_copy = deepcopy(matched_dataset2)
    # get matched RIDs
    matched_subs_dataset1 = (matched_dataset1.RID).values
    matched_subs_dataset2 = (matched_dataset2.RID).values

    # check matched date (visists) all has DX and MMSE
    # otherwise let the date to be nan
    # if all dates of a subject are nan, remove the subject from matched list
    # only forcus on dataset1 is enough
    unqualfied_matched_subs_dataset1 = []
    unqualfied_matched_subs_dataset2 = []
    for i, sub in enumerate(matched_subs_dataset1):
        matched_dates = matched_dataset1.iloc[matched_dataset1[
            matched_dataset1.RID == sub].index, 1:].values[0]
        matched_dates = matched_dates.astype(str)
        matched_dates_val = matched_dates[matched_dates != str(np.nan)]
        for date in matched_dates_val:
            raw_mask = (raw_dataset1.RID == sub) & (
                raw_dataset1.EXAMDATE == date)
            dx = float(raw_dataset1.loc[raw_mask, ['DX']].values[0])
            mmse = float(raw_dataset1.loc[raw_mask, ['MMSE']].values[0])
            # check whether the dx and mmse is nan
            if np.isnan(dx) or np.isnan(mmse):
                # we need to set this date to nan
                # get coordinate
                index = np.where(matched_dates == date)[0]
                matched_dataset1_copy.iloc[matched_dataset1_copy[
                    matched_dataset1_copy.RID == sub].index, index +
                                           1] = np.nan
                matched_dataset2_copy.iloc[matched_dataset1_copy[
                    matched_dataset1_copy.RID == sub].index, index +
                                           1] = np.nan
        # check whether all dates are nan
        matched_dates_copy = matched_dataset1_copy.iloc[matched_dataset1_copy[
            matched_dataset1_copy.RID == sub].index, 1:].values[0]
        matched_dates_copy = matched_dates_copy.astype(str)
        matched_dates_copy = matched_dates_copy[
            matched_dates_copy != str(np.nan)]
        if len(matched_dates_copy) == 0:
            # move this subject to unmatched
            unqualfied_matched_subs_dataset1.append(sub)
            unqualfied_matched_subs_dataset2.append(matched_subs_dataset2[i])
    # get unqualfied indexes
    sorter = np.argsort(matched_subs_dataset1)
    index = sorter[np.searchsorted(
        matched_subs_dataset1, unqualfied_matched_subs_dataset1,
        sorter=sorter)]
    qualfied_matched_subs_dataset1 = np.delete(matched_subs_dataset1, index)
    qualfied_matched_subs_dataset2 = np.delete(matched_subs_dataset2, index)

    # generate unmatched data
    matched_rowIndex_dataset1 = []
    for sub in qualfied_matched_subs_dataset1:
        sub_mask = (raw_dataset1.RID == sub)
        matched_rowIndex_dataset1 += list(raw_dataset1[sub_mask].index.values)
    # drop for dataset1
    unmatched_dataset1 = raw_dataset1.drop(matched_rowIndex_dataset1, axis=0)
    unmatched_dataset1.reset_index(drop=True)
    unmatched_dataset1.to_csv(
        os.path.join(input_path, name1 + '_unmatched.csv'),
        sep=',',
        index=False)

    matched_rowIndex_dataset2 = []
    for sub in qualfied_matched_subs_dataset2:
        sub_mask = (raw_dataset2.RID == sub)
        matched_rowIndex_dataset2 += list(raw_dataset2[sub_mask].index.values)
    # drop for dataset1
    unmatched_dataset2 = raw_dataset2.drop(matched_rowIndex_dataset2, axis=0)
    unmatched_dataset2.reset_index(drop=True)
    unmatched_dataset2.to_csv(
        os.path.join(input_path, name2 + '_unmatched.csv'),
        sep=',',
        index=False)

    # generate matched data
    matched_data_datset1 = pd.DataFrame(columns=cols)
    for sub in qualfied_matched_subs_dataset1:
        matched_dates = matched_dataset1_copy.iloc[matched_dataset1_copy[
            matched_dataset1_copy.RID == sub].index, 1:].values[0]
        matched_dates = matched_dates.astype(str)
        matched_dates = matched_dates[matched_dates != str(np.nan)]
        for date in matched_dates:
            date_mask = (raw_dataset1.RID == sub) & (
                raw_dataset1.EXAMDATE == date)
            row_df = raw_dataset1[date_mask]
            matched_data_datset1 = matched_data_datset1.append(row_df)
    matched_data_datset1 = matched_data_datset1[cols]
    matched_data_datset1.to_csv(
        os.path.join(input_path, name1 + '_matched.csv'), sep=',', index=False)

    matched_data_datset2 = pd.DataFrame(columns=cols)
    for sub in qualfied_matched_subs_dataset2:
        matched_dates = matched_dataset2_copy.iloc[matched_dataset2_copy[
            matched_dataset2_copy.RID == sub].index, 1:].values[0]
        matched_dates = matched_dates.astype(str)
        matched_dates = matched_dates[matched_dates != str(np.nan)]
        for date in matched_dates:
            date_mask = (raw_dataset2.RID == sub) & (
                raw_dataset2.EXAMDATE == date)
            row_df = raw_dataset2[date_mask]
            matched_data_datset2 = matched_data_datset2.append(row_df)
    matched_data_datset2 = matched_data_datset2[cols]
    matched_data_datset2.to_csv(
        os.path.join(input_path, name2 + '_matched.csv'), sep=',', index=False)
