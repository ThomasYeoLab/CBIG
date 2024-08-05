#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import argparse
import numpy as np
import pandas as pd
from statsmodels.multivariate.manova import MANOVA
from evaluation.association_analysis.association_misc import\
    average_rois_multi_folds
from config import global_config
from utils.misc import txt2list, create_folder

import warnings

warnings.filterwarnings('ignore')


def manova_args_parser():
    """
    Feeding parameters
    """
    parser = argparse.ArgumentParser(prog='MANOVArgs')
    parser.add_argument('--input_path', type=str, default='/')
    parser.add_argument('--output_path', type=str, default='/')
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')
    parser.add_argument('--test_name', type=str, default='unmatch2match_test')
    parser.add_argument('--model', type=str, default='TrueBat')
    parser.add_argument('--type', type=str, default='MMSE')
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--norm', action='store_true', default=False)
    parser.add_argument('--demean', action='store_false', default=True)

    args, _ = parser.parse_known_args()
    return args


def nan_manova(df, y_list, x_list):
    """
    Run a MANOVAE analysis using given regressors
    """
    # drop nan value
    df.dropna(subset=y_list, inplace=True)
    model = MANOVA(endog=df[y_list], exog=df[x_list])
    mv_test_results = model.mv_test()
    # extract p values
    model_summary = (mv_test_results.summary_frame['Pr > F']).unstack(level=1)
    # by default use Pillai's trace to compute p values
    p_values = model_summary["Pillai's trace"].values
    assert len(p_values) == len(x_list)

    pillar_traces = (mv_test_results.summary_frame['Value']).unstack(
        level=1)["Pillai's trace"].values
    # convert p values to a dataframe
    p_rows = []
    for i, x in enumerate(x_list):
        row = [x, -1 * np.log10(p_values[i]), pillar_traces[i]]
        p_rows.append(row)
    p_df = pd.DataFrame(data=p_rows, columns=['Beta', 'P', 'PillarTrace'])
    return p_df


def manova_wrapper(args):
    """
    Wrapper function for MANOVA analysis
    """
    ROIs = txt2list(global_config.ROI_features_path)
    gmROIs = txt2list(global_config.gm_ROI_features)

    # Average across 10 folds
    mean_df = average_rois_multi_folds(args.input_path,
                                       [args.test_name + '.csv'],
                                       ROIs,
                                       args.nb_folds,
                                       Demean=args.demean,
                                       ZNorm=args.norm)
    mean_df = mean_df.rename(columns={'EstimatedTotalIntraCranialVol': 'ICV'})
    if args.type == 'MMSE':
        regressors = ['AGE', 'SEX', 'MMSE', 'ICV']
    elif args.type == 'DX':
        regressors = ['AGE', 'SEX', 'MCI', 'AD', 'ICV']
    else:
        raise ValueError('only support MMSE and Diagnosis')

    # run MANOVA
    transROIs = []
    for _, roi in enumerate(gmROIs):
        transROI = roi.replace('-', '')
        transROIs.append(transROI)
        mean_df = mean_df.rename(columns={roi: transROI})

    p_df = nan_manova(mean_df, transROIs, regressors)

    # save results
    if args.norm:
        output_folder = os.path.join(args.output_path,
                                     args.type + '_demean_norm',
                                     args.dataset_pair)
    else:
        output_folder = os.path.join(args.output_path, args.type + '_demean',
                                     args.dataset_pair)
    create_folder(output_folder)
    p_save_name = 'p_' + args.model + '_' + args.dataset_pair + '.csv'
    p_df.to_csv(os.path.join(output_folder, p_save_name), index=False, sep=',')


if __name__ == '__main__':
    manova_wrapper(manova_args_parser())
