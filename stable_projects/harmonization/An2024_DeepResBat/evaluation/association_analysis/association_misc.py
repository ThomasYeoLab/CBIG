#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import copy
import numpy as np
from utils.misc import cat2levels, get_bl_df, csvs2df


def demean_regressors(mean_df):
    """
    Demean regressors
    """
    # perform demean for regressor
    mean_df['AGE'] = (mean_df['AGE'] - mean_df['AGE'].mean())
    mean_df['MMSE'] = (mean_df['MMSE'] - mean_df['MMSE'].mean())
    mean_df['SEX'] = (mean_df['SEX'] - mean_df['SEX'].mean())
    mean_df['EstimatedTotalIntraCranialVol'] = \
        (mean_df['EstimatedTotalIntraCranialVol'] - mean_df[
            'EstimatedTotalIntraCranialVol'].mean())
    # diagnosis
    dx_array = mean_df['DX'].values
    dx_levels = cat2levels(dx_array, isDemean=True)
    mean_df['MCI'] = dx_levels[:, 0]
    mean_df['AD'] = dx_levels[:, 1]

    return mean_df


def average_rois_multi_folds(input_path,
                             file_names,
                             ROIs,
                             nb_folds=10,
                             Demean=True,
                             ZNorm=False):
    """
    Average ROIs from multiple folds
    """
    nb_valid_folds = 0
    first_fold_df = csvs2df(os.path.join(input_path, '0'), file_names)
    first_fold_df = get_bl_df(first_fold_df)
    mean_roi_array = np.zeros((first_fold_df.shape[0], len(ROIs)))
    nb_elements = mean_roi_array.shape[0] * mean_roi_array.shape[1]
    for fold in range(nb_folds):
        fold_input_path = os.path.join(input_path, str(fold))
        full_df = csvs2df(fold_input_path, file_names)
        bl_df = get_bl_df(full_df)
        if fold == 0:
            mean_df = copy.deepcopy(bl_df)
        if np.sum(np.isnan(bl_df[ROIs].values)) < nb_elements:
            mean_roi_array += bl_df[ROIs].values
            nb_valid_folds += 1
        else:
            pass
    mean_df[ROIs] = mean_roi_array / nb_valid_folds
    mean_df['SEX'] -= 1
    if Demean:
        mean_df = demean_regressors(mean_df)
    else:
        dx_array = mean_df['DX'].values
        dx_levels = cat2levels(dx_array, isDemean=False)
        mean_df['MCI'] = dx_levels[:, 0]
        mean_df['AD'] = dx_levels[:, 1]
    if ZNorm:
        rois_mean = mean_df[ROIs].mean()
        rois_std = mean_df[ROIs].std()
        mean_df[ROIs] = (mean_df[ROIs] - rois_mean) / rois_std
    return mean_df


def extract_glm_output(result_df):
    """
    Extract GLM output from statsmodels package
    """
    result = dict()
    keys = list(result_df.index)
    for key in keys:
        result[key] = dict()
        coef = result_df.loc[key, 'coef']
        std_err = result_df.loc[key, 'std err']
        z = result_df.loc[key, 'z']
        p = result_df.loc[key, 'P>|z|']
        result[key]['coef'] = coef
        result[key]['std_err'] = std_err
        result[key]['z'] = z
        result[key]['p'] = p
    return result
