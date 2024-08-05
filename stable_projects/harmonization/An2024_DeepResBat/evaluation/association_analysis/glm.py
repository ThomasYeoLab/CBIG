#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import argparse
import pandas as pd
import statsmodels.formula.api as smf
from evaluation.association_analysis.association_misc import\
    extract_glm_output, average_rois_multi_folds
from config import global_config
from utils.misc import txt2list, create_folder

formula_dict = {
    'MMSE': " ~ AGE + SEX + MMSE + ICV",
    'DX': " ~ AGE + SEX + MCI + AD + ICV",
    "LeftHipp": " ~ AGE + SEX + MMSE + ICV"
}


def nan_glm(df, roi, glm_formula):
    """
    Fit a GLM model with given formula
    """
    # drop nan value
    df.dropna(subset=[roi], inplace=True)
    model = smf.glm(formula=glm_formula, data=df).fit()
    # extract results
    model_summary = model.summary()
    results_as_html = model_summary.tables[1].as_html()
    result_df = pd.read_html(results_as_html, header=0, index_col=0)[0]
    result = extract_glm_output(result_df)

    return result


def glm_runner(df, gmROIs, formula, regresssors):
    """
    Run GLM for all grey matter ROIs
    """
    beta_rows, z_rows, p_rows = [], [], []
    for roi in gmROIs:
        transROI = roi.replace('-', '')
        df = df.rename(columns={roi: transROI})
        glm_formula = transROI + formula
        glm_results = nan_glm(df, transROI, glm_formula)
        # beta
        beta_row = [roi] + [
            glm_results[regressor]['coef'] for regressor in regresssors
        ]
        beta_rows.append(beta_row)
        # z statisctics
        z_row = [roi] + [
            glm_results[regressor]['z'] for regressor in regresssors
        ]
        z_rows.append(z_row)
        # p value
        p_row = [roi] + [
            glm_results[regressor]['p'] for regressor in regresssors
        ]
        p_rows.append(p_row)
    cols = ['ROI'] + regresssors
    beta_df = pd.DataFrame(data=beta_rows, columns=cols)
    z_df = pd.DataFrame(data=z_rows, columns=cols)
    p_df = pd.DataFrame(data=p_rows, columns=cols)
    return beta_df, z_df, p_df


def glm_args_parser():
    """
    Feeding parameters
    """
    parser = argparse.ArgumentParser(prog='GLMArgs')
    parser.add_argument('--input_path', type=str, default='/')
    parser.add_argument('--output_path', type=str, default='/')
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')
    parser.add_argument('--test_name', type=str, default='unmatch2match_test')
    parser.add_argument('--model', type=str, default='TrueBat')
    parser.add_argument('--type', type=str, default='MMSE')
    parser.add_argument('--norm', action='store_true', default=False)
    parser.add_argument('--demean', action='store_false', default=True)

    args, _ = parser.parse_known_args()
    return args


def glm_wrapper(args):
    """
    Wrapper function for performing GLM analysis
    """
    ROIs = txt2list(global_config.ROI_features_path)
    if args.type == 'LeftHipp':
        gmROIs = ['Left-Hippocampus']
    else:
        gmROIs = txt2list(global_config.gm_ROI_features)
    # Average test data across 10 folds
    mean_df = average_rois_multi_folds(args.input_path,
                                       [args.test_name + '.csv'],
                                       ROIs,
                                       Demean=args.demean,
                                       ZNorm=args.norm)
    mean_df = mean_df.rename(columns={'EstimatedTotalIntraCranialVol': 'ICV'})
    if args.type == 'MMSE':
        regressors = ['Intercept', 'AGE', 'SEX', 'ICV', 'MMSE']
    elif args.type == 'DX':
        regressors = ['Intercept', 'AGE', 'SEX', 'ICV', 'MCI', 'AD']
    elif args.type == 'LeftHipp':
        regressors = ['Intercept', 'AGE', 'SEX', 'ICV', 'MMSE']
    else:
        raise ValueError('only support MMSE and Diagnosis')
    # GLM
    beta_df, z_df, p_df = glm_runner(mean_df, gmROIs, formula_dict[args.type],
                                     regressors)
    if args.demean:
        if args.norm:
            output_folder = os.path.join(args.output_path,
                                         args.type + '_demean_norm',
                                         args.dataset_pair)
            create_folder(output_folder)
        else:
            output_folder = os.path.join(args.output_path,
                                         args.type + '_demean',
                                         args.dataset_pair)
            create_folder(output_folder)
    else:
        raise ValueError('only support demean for regressors')
    # save glm output
    beta_save_name = 'beta_' + args.model + '_' + args.dataset_pair + '.csv'
    beta_df.to_csv(os.path.join(output_folder, beta_save_name),
                   index=False,
                   sep=',')
    z_save_name = 'z_' + args.model + '_' + args.dataset_pair + '.csv'
    z_df.to_csv(os.path.join(output_folder, z_save_name), index=False, sep=',')
    p_save_name = 'p_' + args.model + '_' + args.dataset_pair + '.csv'
    p_df.to_csv(os.path.join(output_folder, p_save_name), index=False, sep=',')


if __name__ == '__main__':
    glm_wrapper(glm_args_parser())
