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
from utils.misc import create_folder, load_pkl, save_pkl


def pick_well_matched_pairs(args,
                            step4_output_path,
                            step5_output_path,
                            name1='MACC',
                            name2='ADNI'):
    """
    Pick time points with well matched MMSE

    Args:
        args (tuple): Parameters
        step4_output_path (str): Path for step4 output
        step5_output_path (str): Path for step5 output
        name1 (str, optional): Name for dataset1. Defaults to 'MACC'.
        name2 (str, optional): Name for dataset2. Defaults to 'ADNI'.
    """
    step4_threshold_path = os.path.join(step4_output_path,
                                        'Threshold_' + str(args.threshold))
    step5_threshold_path = os.path.join(step5_output_path,
                                        'Threshold_' + str(args.threshold))
    create_folder(step5_threshold_path)
    # read matched_df
    matched_df = pd.read_csv(os.path.join(step4_threshold_path, 'matched.csv'))
    raw_dataset1 = pd.read_csv(os.path.join(args.MACC_data_path))
    raw_dataset2 = pd.read_csv(os.path.join(args.ADNI_data_path))
    subs_dataset1 = np.unique(matched_df.MACC_RID)
    cols_dataset1 = [
        'MACC_1', 'MACC_2', 'MACC_3', 'MACC_4', 'MACC_5', 'MACC_6'
    ]
    cols_dataset1 = cols_dataset1[:args.threshold]
    cols_dataset2 = [
        'ADNI_1', 'ADNI_2', 'ADNI_3', 'ADNI_4', 'ADNI_5', 'ADNI_6'
    ]
    cols_dataset2 = cols_dataset2[:args.threshold]

    # pick well matched time points
    for sub_dataset1 in subs_dataset1:
        # get matched subject pair
        matched_row = matched_df[(matched_df.MACC_RID == sub_dataset1)]
        dates_dataset1 = matched_row[cols_dataset1].values.astype(str)
        dates_dataset1 = dates_dataset1[dates_dataset1 != str(np.nan)]

        sub_dataset2 = matched_row['ADNI_RID'].values[0]
        dates_dataset2 = matched_row[cols_dataset2].values.astype(str)
        dates_dataset2 = dates_dataset2[dates_dataset2 != str(np.nan)]

        assert len(dates_dataset1) == len(
            dates_dataset2), 'Ineqaul matched TPs'

        # Test whether MMSE differnce is larger than args.epslion
        nb_visits = len(dates_dataset1)
        for i in range(nb_visits):
            date_dataset1 = dates_dataset1[i]
            mask_dataset1 = (raw_dataset1.RID == sub_dataset1) & (
                raw_dataset1.EXAMDATE == date_dataset1)
            if raw_dataset1[mask_dataset1]['MMSE'].values.shape[0] >= 1:
                mmse_dataset1 = raw_dataset1[mask_dataset1]['MMSE'].values[0]
            else:
                mmse_dataset1 = np.nan

            date_dataset2 = dates_dataset2[i]
            mask_dataset2 = (raw_dataset2.RID == sub_dataset2) & (
                raw_dataset2.EXAMDATE == date_dataset2)
            if raw_dataset2[mask_dataset2]['MMSE'].values.shape[0] >= 1:
                mmse_dataset2 = raw_dataset2[mask_dataset2]['MMSE'].values[0]
            else:
                mmse_dataset2 = np.nan

            # checkwether MMSE difference is within epslion
            if not np.isnan(mmse_dataset1) and not np.isnan(mmse_dataset2):
                if abs(mmse_dataset1 - mmse_dataset2) > args.mmse_eps:
                    # replace matched date with NaN
                    matched_df.iloc[matched_df[(
                        matched_df.MACC_RID == sub_dataset1)].index, i +
                                    1] = np.nan
                    matched_df.iloc[matched_df[(
                        matched_df.MACC_RID == sub_dataset1)].index, i + 2 +
                                    args.threshold] = np.nan
            else:
                matched_df.iloc[matched_df[(matched_df.MACC_RID == sub_dataset1
                                            )].index, i + 1] = np.nan
                matched_df.iloc[matched_df[(
                    matched_df.MACC_RID == sub_dataset1)].index, i + 2 +
                                args.threshold] = np.nan
    # we need to clean matched_df
    matched_df = clean_dropped_matched(matched_df, args.threshold)
    matched_df = rm_empty_rows(matched_df, cols_dataset1)
    matched_df = rm_empty_cols(matched_df, cols_dataset1, cols_dataset2)
    # save the matched_df
    matched_df.to_csv(
        os.path.join(step5_threshold_path, 'picked_matched.csv'),
        sep=',',
        index=False)


def clean_dropped_matched(matched_df, threshold):
    """
    Clean dropped matched df by moving date to deleted date

    Args:
        matched_df (class DataFrame): Matched dataframe
        threshold (int): Threshold to beigin matching
    """
    for i in range(len(matched_df)):
        for k in range(threshold):
            macc_date = matched_df.iloc[i, k + 1]
            if isinstance(macc_date, float):
                for j in range(k + 1, threshold):
                    if isinstance(matched_df.iloc[i, j + 1], str):
                        matched_df.iloc[i, k + 1] = matched_df.iloc[i, j + 1]
                        matched_df.iloc[i, j + 1] = np.nan
                        matched_df.iloc[i, k + 2 + threshold] = \
                            matched_df.iloc[i, j + 2 + threshold]
                        matched_df.iloc[i, j + 2 + threshold] = np.nan
                        continue
    return matched_df


def rm_empty_rows(matched_df, cols_dataset1):
    """
    Remove rows where all timepoints are poorly matched

    Args:
        matched_df (class DataFrame): Matched dataframe
        cols_dataset1 (list): Columns for dataset1
    """
    _matched_df = deepcopy(matched_df)
    for i in range(len(matched_df)):
        row = matched_df.iloc[i]
        dates = row[cols_dataset1].values
        dates = dates.astype(str)
        dates = dates[dates != str(np.nan)]
        if dates.shape[0] == 0:
            # drop!
            _matched_df.drop(i, inplace=True)

    return _matched_df


def rm_empty_cols(matched_df, cols_dataset1, cols_dataset2):
    """
    Remove columns where all timepoints are poorly matched

    Args:
        matched_df (class DataFrame): Matched dataframe
        cols_dataset1 (list): Columns for dataset1
        cols_dataset2 (list): Columns for dataset1
    """
    _matched_df = deepcopy(matched_df)
    empty_cols = []
    for i in range(len(cols_dataset1)):
        series = matched_df[cols_dataset1[i]].values
        series = series.astype(str)
        series = series[series != str(np.nan)]
        if series.shape[0] == 0:
            empty_cols.append(cols_dataset1[i])
            empty_cols.append(cols_dataset2[i])

    _matched_df.drop(empty_cols, axis=1, inplace=True)

    return _matched_df


def split_matched_csv(args, bin, step5_output_path, name1='MACC',
                      name2='ADNI'):
    """
    Split the matched csv file into MACC(AIBL) and ADNI

    Args:
        args (tuple): Parameters
        bin (int): Bin
        step5_output_path (str): Path for step4 output
        name1 (str, optional): Name for dataset1. Defaults to 'MACC'.
        name2 (str, optional): Name for dataset2. Defaults to 'ADNI'.
    """
    output_path = os.path.join(step5_output_path,
                               'Threshold_' + str(args.threshold))
    # read matched csv path
    matched_df = pd.read_csv(os.path.join(output_path, 'picked_matched.csv'))
    nb_tps = int(matched_df.shape[1] / 2 - 1)
    # MACC (AIBL)
    cols_dataset1 = []
    cols_dataset1.append('MACC_RID')
    for t in range(1, nb_tps + 1):
        tp_name = name1 + '_' + str(t)
        cols_dataset1.append(tp_name)
    matched_dataset1 = pd.read_csv(
        os.path.join(output_path, 'picked_matched.csv'), usecols=cols_dataset1)
    new_cols_dataset1 = []
    new_cols_dataset1.append('RID')
    for t in range(1, nb_tps + 1):
        tp_name = str(t)
        new_cols_dataset1.append(tp_name)
    matched_dataset1.columns = new_cols_dataset1
    # ADNI
    cols_dataset2 = []
    cols_dataset2.append('ADNI_RID')
    for t in range(1, nb_tps + 1):
        tp_name = name2 + '_' + str(t)
        cols_dataset2.append(tp_name)

    matched_dataset2 = pd.read_csv(
        os.path.join(output_path, 'picked_matched.csv'), usecols=cols_dataset2)
    new_cols_dataset2 = []
    new_cols_dataset2.append('RID')
    for t in range(1, nb_tps + 1):
        tp_name = str(t)
        new_cols_dataset2.append(tp_name)

    matched_dataset2.columns = new_cols_dataset2

    # we need to control #matched subjects
    raw_dataset1 = pd.read_csv(
        os.path.join(args.checkpoint_path, args.matching_pair,
                     'matching_' + str(args.nb_bins) + 'BINs',
                     'BIN_' + str(bin), 'MACC_' + str(bin) + '_bin.csv'))
    raw_dataset2 = pd.read_csv(args.ADNI_data_path)
    nb_bin_subjects = len(np.unique(raw_dataset1.RID))
    if matched_dataset1.shape[0] <= round(args.match_ratio * nb_bin_subjects):
        matched_dataset1.to_csv(
            os.path.join(output_path, 'picked_' + name1 + '.csv'),
            sep=',',
            index=False)

        matched_dataset2.to_csv(
            os.path.join(output_path, 'picked_' + name2 + '.csv'),
            sep=',',
            index=False)
        matched_subs_dataset1 = np.unique(matched_dataset1.RID)
        matched_subs_dataset2 = np.unique(matched_dataset2.RID)
    else:
        # we need to drop some subjects
        matched_cost = gen_matched_subjects_cost_table(
            args, raw_dataset1, raw_dataset2, matched_dataset1,
            matched_dataset2, output_path)
        rm_macc_subs = matched_cost.iloc[int(args.match_ratio *
                                             nb_bin_subjects):, 0].values
        rm_adni_subs = matched_cost.iloc[int(args.match_ratio *
                                             nb_bin_subjects):, 1].values
        rm_high_cost_pairs(
            rm_macc_subs, matched_dataset1,
            os.path.join(output_path, 'picked_' + name1 + '.csv'))
        rm_high_cost_pairs(
            rm_adni_subs, matched_dataset2,
            os.path.join(output_path, 'picked_' + name2 + '.csv'))
        dataset1 = pd.read_csv(
            os.path.join(output_path, 'picked_' + name1 + '.csv'))
        dataset2 = pd.read_csv(
            os.path.join(output_path, 'picked_' + name2 + '.csv'))
        matched_subs_dataset1 = np.unique(dataset1.RID)
        matched_subs_dataset2 = np.unique(dataset2.RID)

    return matched_subs_dataset1, matched_subs_dataset2


def gen_matched_subjects_cost_table(args, raw_dataset1, raw_dataset2,
                                    picked_MACC, picked_ADNI, output_path):
    """
    Generate cost table for matched subejcts to drop high cost pairs

    Args:
        args (tuple): Parameters
        raw_dataset1 (class DataFrame): Dataframe for raw dataset1
        raw_dataset2 (class DataFrame): Dataframe for raw dataset1
        picked_MACC (class DataFrame): Dataframe for picked dataset1
        picked_ADNI (class DataFrame): Dataframe for picked dataset2
        output_path (str): Path for saving output
    """
    matched_macc_subjects = np.unique(picked_MACC.RID)
    cost_rows = []
    for macc_sub in matched_macc_subjects:
        cost_row = []
        # get matched_macc dates and matched_adni sub and dates
        sub_mask = (picked_MACC.RID == macc_sub)
        macc_dates = picked_MACC.iloc[picked_MACC[sub_mask].
                                      index, 1:].values[0]
        macc_dates = macc_dates.astype(str)
        macc_dates = macc_dates[macc_dates != str(np.nan)]
        adni_sub = picked_ADNI.iloc[picked_MACC[sub_mask].index, 0].values[0]
        adni_dates = picked_ADNI.iloc[picked_MACC[sub_mask].
                                      index, 1:].values[0]
        adni_dates = adni_dates.astype(str)
        adni_dates = adni_dates[adni_dates != str(np.nan)]
        assert len(adni_dates) == len(adni_dates), 'Wrong matching case'

        cost, mmse_cost = matching_cost(args, raw_dataset1, raw_dataset2,
                                        macc_sub, macc_dates, adni_sub,
                                        adni_dates)
        cost_row.append(macc_sub)
        cost_row.append(adni_sub)
        cost_row.append(cost)
        cost_row.append(mmse_cost)
        cost_rows.append(cost_row)

    # generate dataframe
    matched_cost = pd.DataFrame(
        data=cost_rows, columns=['MACC_RID', 'ADNI_RID', 'Cost', 'MMSE_cost'])
    # sort according to cost
    matched_cost.sort_values(
        by=['Cost'], inplace=True, ascending=True, ignore_index=True)
    matched_cost.to_csv(
        os.path.join(output_path, 'cost.csv'), sep=',', index=False)

    return matched_cost


def matching_cost(args, raw_dataset1, raw_dataset2, macc_sub, macc_dates,
                  adni_sub, adni_dates):
    """
    Calculate matching cost for matched subject pairs

    Args:
        args (tuple): Parameters
        raw_dataset1 (class DataFrame): Dataframe for raw dataset1
        raw_dataset2 (class DataFrame): Dataframe for raw dataset1
        macc_sub (list): List of MACC subject RIDs
        macc_dates (list): List of MACC dates
        adni_sub (list): List of ADNI subject RIDs
        adni_dates (list): List of ADNI dates
    """
    cost = 0
    mmse_cost = 0
    nb_dates = len(macc_dates)
    for i, macc_date in enumerate(macc_dates):
        macc_mask = (raw_dataset1.RID == macc_sub) & (
            raw_dataset1.EXAMDATE == macc_date)
        adni_mask = (raw_dataset2.RID == adni_sub) & (
            raw_dataset2.EXAMDATE == adni_dates[i])
        # cost for age
        macc_age = raw_dataset1.loc[macc_mask, ['AGE']].values[0][0]
        adni_age = raw_dataset2.loc[adni_mask, ['AGE']].values[0][0]
        cost += matching_cost_one_measure(args.age_penalty, args.NANpenalty,
                                          macc_age, adni_age)
        # cost for sex
        macc_sex = raw_dataset1.loc[macc_mask, ['SEX']].values[0][0]
        adni_sex = raw_dataset2.loc[adni_mask, ['SEX']].values[0][0]
        cost += matching_cost_one_measure(args.sex_penalty, args.NANpenalty,
                                          macc_sex, adni_sex)
        # cost for dx
        macc_dx = raw_dataset1.loc[macc_mask, ['DX']].values[0][0]
        adni_dx = raw_dataset2.loc[adni_mask, ['DX']].values[0][0]
        cost += matching_cost_one_measure(args.dx_penalty, args.NANpenalty,
                                          macc_dx, adni_dx)
        # cost for mmse
        macc_mmse = raw_dataset1.loc[macc_mask, ['MMSE']].values[0][0]
        adni_mmse = raw_dataset2.loc[adni_mask, ['MMSE']].values[0][0]
        mmse_cost += matching_cost_one_measure(
            args.mmse_penalty, args.NANpenalty, macc_mmse, adni_mmse)
        cost += mmse_cost
    return cost / nb_dates, mmse_cost / nb_dates


def matching_cost_one_measure(penalty, NANpenalty, macc, adni):
    """
    Calulate matching cost for one measure for matched time points

    Args:
        penalty (int): Penalty for punishing badly matched
        NANpenalty (int): Penalty for punishing missing data
        macc (ndarray): Vector for MACC
        adni (ndarray): Vector for ADNI
    """
    if np.isnan(macc):
        return NANpenalty
    elif np.isnan(adni):
        return NANpenalty
    else:
        return penalty * np.abs(adni - macc)


def rm_high_cost_pairs(rm_subs, picked_df, save_path):
    """
    Remove matched pairs with relatively high cost

    Args:
        rm_subs (list): List of subjectw with high costs
        picked_df (class DataFrame): Picked matched dataframe
        save_path (str): Path for saving output
    """
    row_indexes = []
    for sub in rm_subs:
        index = picked_df[picked_df.RID == sub].index.values[0]
        row_indexes.append(index)
    # drop
    controled_df = picked_df.drop(row_indexes, axis=0)
    controled_df.reset_index(drop=True, inplace=True)
    # note that we need to check whether having empty cols
    cols = list(controled_df.columns)
    df = deepcopy(controled_df)
    empty_cols = []
    for i in range(len(cols)):
        series = controled_df[cols[i]].values
        series = series.astype(str)
        series = series[series != str(np.nan)]
        if series.shape[0] == 0:
            empty_cols.append(cols[i])
    df.drop(empty_cols, axis=1, inplace=True)
    # save
    df.to_csv(save_path, sep=',', index=False)


def update_next_round_data(args,
                           bin,
                           matched_subs_dataset1,
                           matched_subs_dataset2,
                           name1='MACC',
                           name2='ADNI'):
    """
    Remove the matched MACC subjects and ADNI subjects in next round

    Args:
        args (tuple): Parameters
        bin (int): Bin
        matched_subs_dataset1 (list): List of matched subjects RIDs of dataset1
        matched_subs_dataset2 (list): List of matched subjects RIDs of dataset1
        name1 (str, optional): Name for dataset1. Defaults to 'MACC'.
        name2 (str, optional): Name for dataset2. Defaults to 'ADNI'.
    """
    curr_round_data_path = os.path.join(
        args.checkpoint_path, args.matching_pair,
        'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
        'round_' + str(args.round))
    next_round_data_path = os.path.join(
        args.checkpoint_path, args.matching_pair,
        'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
        'round_' + str(args.round + 1))
    create_folder(next_round_data_path)
    save_name1 = name1 + '_' + 'AGE' + str(args.age_penalty) + '_' + \
        'SEX' + str(args.sex_penalty) + '_' + \
        'DX' + str(args.dx_penalty) + '_' \
        'MMSE' + str(args.mmse_penalty) + '.pkl'
    save_name2 = name2 + '_' + 'AGE' + str(args.age_penalty) + '_' + \
        'SEX' + str(args.sex_penalty) + '_' + \
        'DX' + str(args.dx_penalty) + '_' \
        'MMSE' + str(args.mmse_penalty) + '.pkl'
    # load data
    if args.round == 1 and bin == 0:
        comb_dataset1 = load_pkl(
            os.path.join(curr_round_data_path, name1 + '.pkl'))
        comb_dataset2 = load_pkl(
            os.path.join(curr_round_data_path, name2 + '.pkl'))
    elif args.round == 1 and bin > 0:
        comb_dataset1 = load_pkl(
            os.path.join(curr_round_data_path, name1 + '.pkl'))
        comb_dataset2 = load_pkl(
            os.path.join(curr_round_data_path, save_name2))
    else:
        comb_dataset1 = load_pkl(
            os.path.join(curr_round_data_path, save_name1))
        comb_dataset2 = load_pkl(
            os.path.join(curr_round_data_path, save_name2))

    # deep copy
    _comb_dataset1 = deepcopy(comb_dataset1)
    _comb_dataset2 = deepcopy(comb_dataset2)

    for sub_dataset1 in _comb_dataset1.keys():
        if sub_dataset1 in matched_subs_dataset1:
            del comb_dataset1[sub_dataset1]

    for sub_dataset2 in _comb_dataset2.keys():
        if sub_dataset2 in matched_subs_dataset2:
            del comb_dataset2[sub_dataset2]

    # save updated within dataset combs
    save_pkl(comb_dataset1, os.path.join(next_round_data_path, save_name1))
    save_pkl(comb_dataset2, os.path.join(next_round_data_path, save_name2))
