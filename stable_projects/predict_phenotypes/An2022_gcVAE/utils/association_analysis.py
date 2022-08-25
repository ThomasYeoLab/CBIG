#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import numpy as np
import pandas as pd
import scipy.stats
from config import global_config
from utils.misc import txt2list, csvs2df


def nan_corr(brain_volume, variable):
    """
    Calculate the correlation between a brain volume and a variable

    Args:
        brain_volume (ndarray): Vector for <ROI> volume
        variable (ndarray): Vector for interested variable
    """
    assert brain_volume.shape[0] == variable.shape[0], 'Not equal length!'
    assert brain_volume.shape[1] == 1, 'brain volume should be N by 1'
    assert variable.shape[1] == 1, 'variable should be N by 1'

    # remove nan
    nan_index_brain_volume = np.isnan(brain_volume)
    nan_index_variable = np.isnan(variable)
    nan_index = np.logical_or(nan_index_brain_volume, nan_index_variable)
    nonan_brain_volume = brain_volume[~nan_index]
    nonan_variable = variable[~nan_index]

    # compute correlation
    r, _ = scipy.stats.pearsonr(nonan_brain_volume, nonan_variable)

    return r


def nan_etasq(brain_volume, variable, nb_levels):
    """
    Calculate the eta^2 between a brain volume and a variable

    Args:
        brain_volume (ndarray): Vector for <ROI> volume
        variable (ndarray): Vector for interested variable
        nb_levels (int): Levels for categorarical variable
    """
    assert brain_volume.shape[0] == variable.shape[0], 'Not equal length!'
    assert brain_volume.shape[1] == 1, 'brain volume should be N by 1'
    assert variable.shape[1] == 1, 'variable should be N by 1'

    # remove nan
    nan_index_brain_volume = np.isnan(brain_volume)
    nan_index_variable = np.isnan(variable)
    nan_index = np.logical_or(nan_index_brain_volume, nan_index_variable)
    nonan_brain_volume = brain_volume[~nan_index]
    nonan_variable = variable[~nan_index]
    # compute grand_mean
    grand_mean = np.mean(nonan_brain_volume)
    # compute group_means
    group_means = dict()
    for i in range(nb_levels):
        group_means[i] = np.mean(nonan_brain_volume[nonan_variable == i])
    # compute between-groups sum of squares SSB
    # print(grand_mean, group_means)
    SSB = 0
    for i in range(nb_levels):
        N_i = nonan_brain_volume[nonan_variable == i].shape[0]
        sq = (group_means[i] - grand_mean)**2
        SSB += N_i * sq
    # compute within-groups sum of square SSW
    SSW = 0
    for i in range(nb_levels):
        nonan_brain_volume_i = nonan_brain_volume[nonan_variable == i]
        N_i = nonan_brain_volume_i.shape[0]
        for j in range(N_i):
            sq = (nonan_brain_volume_i[j] - group_means[i])**2
            SSW += sq
    SST = SSB + SSW

    return SSB / SST


def corr_1df(df, cov, ROIs):
    """
    Compute correlation between cov and ROIs for a DataFrame

    Args:
        df (class DataFrame): Brain ROI volumes
        cov (ndarray): Vector for cov
        ROIs (list): List of ROIs
    """
    rows = []
    nb_visits = df.shape[0]
    assert cov in ['AGE', 'MMSE'], "cov not in ['AGE','MMSE']"
    for roi in ROIs:
        roi_array = (df[roi].values).reshape((nb_visits, 1))
        cov_array = (df[cov].values).reshape((nb_visits, 1))
        r = nan_corr(roi_array, cov_array)
        row = [roi, r]
        rows.append(row)
    r_df_cols = ['ROI', cov]
    r_df = pd.DataFrame(data=rows, columns=r_df_cols)
    return r_df


def etasq_1df(df, cov, ROIs):
    """
    Compute eta^2 between cov and ROIs for a DataFrame

    Args:
        df (class DataFrame): Brain ROI volumes
        cov (ndarray): Vector for cov
        ROIs (list): List of ROIs
    """
    rows = []
    nb_visits = df.shape[0]
    assert cov in ['SEX', 'DX'], "cov not in ['SEX','DX']"
    if cov == 'SEX':
        nb_levels = 2
        df.loc[:, ['SEX']] -= 1
    else:
        nb_levels = 3
    for roi in ROIs:
        roi_array = (df[roi].values).reshape((nb_visits, 1))
        cov_array = (df[cov].values).reshape((nb_visits, 1))
        eta_sq = nan_etasq(roi_array, cov_array, nb_levels)
        row = [roi, eta_sq]
        rows.append(row)
    etasq_df_cols = ['ROI', cov]
    etasq_df = pd.DataFrame(data=rows, columns=etasq_df_cols)
    return etasq_df


def corr_10folds(input_path, output_path, ROIs, model, dataset_pair):
    """
    Average correlation across 10 folds

    Args:
        input_path (str): Path for input data
        output_path (str): Pathf for output data
        ROIs (list): List of ROIs
        model (str): Harmonization model name
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
    """
    nb_folds = 10
    nb_rois = 87  # grey matter ROIs
    COVs = ['AGE', 'MMSE']
    mean_r_array = np.zeros((nb_rois, len(COVs)))
    models = ['Unharm', 'ComBat', 'cVAE', 'gcVAE']
    assert model in models, "model not in ['Unharm','ComBat','cVAE','gcVAE']"
    # average harmonized ROIs across 10 folds
    for fold in range(nb_folds):
        fold_input_path = os.path.join(input_path, str(fold))
        if model == 'Unharm':
            csvs = ['unmatch2match_test.csv']
        elif model == 'ComBat':
            csvs = ['unmatch2match_test.csv']
        else:
            csvs = ['unmatch2match_test-intermediate.csv']
        df = csvs2df(fold_input_path, csvs)
        if fold == 0:
            mean_roi_array = np.zeros((df.shape[0], nb_rois))
        mean_roi_array += df[ROIs].values / nb_folds
    mean_roi_df = pd.DataFrame(data=mean_roi_array, columns=ROIs)
    mean_roi_df['AGE'] = df['AGE'].values
    mean_roi_df['MMSE'] = df['MMSE'].values
    age_fold_r_df = corr_1df(mean_roi_df, COVs[0], ROIs)
    mmse_fold_r_df = corr_1df(mean_roi_df, COVs[1], ROIs)
    fold_r_df = pd.concat([age_fold_r_df, mmse_fold_r_df.iloc[:, 1:]], axis=1)
    mean_r_array = fold_r_df[COVs].values
    mean_r_df_cols = ['ROI'] + COVs
    mean_r_df = pd.DataFrame(columns=mean_r_df_cols)
    mean_r_df.iloc[:, 0] = ROIs
    mean_r_df.iloc[:, 1:] = mean_r_array
    mean_r_df_save_name = dataset_pair + '_' + model + '_mean_r_AGE_MMSE.csv'
    mean_r_df.to_csv(
        os.path.join(output_path, mean_r_df_save_name), index=False)


def etasq_10folds(input_path, output_path, ROIs, model, dataset_pair):
    """
    Average eta^2 across 10 folds

    Args:
        input_path (str): Path for input data
        output_path (str): Pathf for output data
        ROIs (list): List of ROIs
        model (str): Harmonization model name
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
    """
    nb_folds = 10
    nb_rois = 87  # grey matter ROIs
    COVs = ['SEX', 'DX']
    mean_etqsq_array = np.zeros((nb_rois, len(COVs)))
    models = ['Unharm', 'ComBat', 'cVAE', 'gcVAE']
    assert model in models, "model not in ['Unharm','ComBat','cVAE','gcVAE']"
    # average across 10 folds
    for fold in range(nb_folds):
        fold_input_path = os.path.join(input_path, str(fold))
        if model == 'Unharm':
            csvs = ['unmatch2match_test.csv']
        elif model == 'ComBat':
            csvs = ['unmatch2match_test.csv']
        else:
            csvs = ['unmatch2match_test-intermediate.csv']
        df = csvs2df(fold_input_path, csvs)
        if fold == 0:
            mean_roi_array = np.zeros((df.shape[0], nb_rois))
        mean_roi_array += df[ROIs].values / nb_folds
    mean_roi_df = pd.DataFrame(data=mean_roi_array, columns=ROIs)
    mean_roi_df['SEX'] = df['SEX'].values
    mean_roi_df['DX'] = df['DX'].values
    sex_fold_etasq_df = etasq_1df(mean_roi_df, COVs[0], ROIs)
    dx_fold_etasq_df = etasq_1df(mean_roi_df, COVs[1], ROIs)
    fold_etasq_df = pd.concat(
        [sex_fold_etasq_df, dx_fold_etasq_df.iloc[:, 1:]], axis=1)
    mean_etqsq_array += fold_etasq_df[COVs].values
    mean_etasq_df_cols = ['ROI'] + COVs
    mean_etasq_df = pd.DataFrame(columns=mean_etasq_df_cols)
    mean_etasq_df.iloc[:, 0] = ROIs
    mean_etasq_df.iloc[:, 1:] = mean_etqsq_array
    mean_etasq_df_save_name =\
        dataset_pair + '_' + model + '_mean_etasq_SEX_DX.csv'
    mean_etasq_df.to_csv(
        os.path.join(output_path, mean_etasq_df_save_name), index=False)


def main():
    """
    Main function for performing association analysis
    """
    # get root path
    root_path = os.path.abspath(
        os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
    data_path = os.path.join(root_path, 'data')
    results_path = os.path.join(root_path, 'results', 'unmatch2match',
                                'association_analysis')
    ROIs = txt2list(global_config.gm_ROI_features)
    # for ADNI-AIBL
    dataset_pairs = ['ADNI-AIBL', 'ADNI-MACC']
    for dataset_pair in dataset_pairs:
        # for Unharm
        corr_10folds(
            os.path.join(data_path, 'splits', dataset_pair),
            os.path.join(results_path, dataset_pair), ROIs, 'Unharm',
            dataset_pair)
        etasq_10folds(
            os.path.join(data_path, 'splits', dataset_pair),
            os.path.join(results_path, dataset_pair), ROIs, 'Unharm',
            dataset_pair)
        # for ComBat harmonization
        corr_10folds(
            os.path.join(data_path, 'unmatch2match', 'harm_output',
                         dataset_pair, 'ComBat', 'AGE_SEX_noref'),
            os.path.join(results_path, dataset_pair), ROIs, 'ComBat',
            dataset_pair)
        etasq_10folds(
            os.path.join(data_path, 'unmatch2match', 'harm_output',
                         dataset_pair, 'ComBat', 'AGE_SEX_noref'),
            os.path.join(results_path, dataset_pair), ROIs, 'ComBat',
            dataset_pair)
        # for cVAE harmonization
        corr_10folds(
            os.path.join(data_path, 'unmatch2match', 'harm_output',
                         dataset_pair, 'cVAE'),
            os.path.join(results_path, dataset_pair), ROIs, 'cVAE',
            dataset_pair)
        etasq_10folds(
            os.path.join(data_path, 'unmatch2match', 'harm_output',
                         dataset_pair, 'cVAE'),
            os.path.join(results_path, dataset_pair), ROIs, 'cVAE',
            dataset_pair)
        # for gcVAE harmonization
        corr_10folds(
            os.path.join(data_path, 'unmatch2match', 'harm_output',
                         dataset_pair, 'gcVAE'),
            os.path.join(results_path, dataset_pair), ROIs, 'gcVAE',
            dataset_pair)
        etasq_10folds(
            os.path.join(data_path, 'unmatch2match', 'harm_output',
                         dataset_pair, 'gcVAE'),
            os.path.join(results_path, dataset_pair), ROIs, 'gcVAE',
            dataset_pair)


if __name__ == '__main__':
    main()
