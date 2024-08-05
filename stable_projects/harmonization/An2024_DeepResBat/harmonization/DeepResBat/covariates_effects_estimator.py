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
from typing import List
from datetime import datetime
from config import global_config
from model.XGB import XGBWrapper
from harmonization.DeepResBat.harmonization_misc import \
    reliability_voting, global_checker, save_covs_effects, regress_out_site,\
    linear_covariates_effects_estimator, data_proc_estimate_covar_effects
from utils.misc import\
    txt2list, create_folder, save_pkl, load_pkl, list2txt


covs_ref_dict = {
    "AGE": "regress",
    "SEX": "binary",
    "MMSE": "regress",
    "DX": "multi"
}


def covariates_effects_estimator_parser():
    """
    Parameters for non-lienar predictor F(X)
    """
    parser = argparse.ArgumentParser(prog='CovPredArgs')
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--checkpoint_path', type=str, default='/')
    parser.add_argument('--hyper_params_path', type=str, default='/')
    parser.add_argument('--output_path', type=str, default='/')
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--model', type=str, default='DeepResBat')
    parser.add_argument('--sufix', type=str, default='_F')
    parser.add_argument('--train_name',
                        type=str,
                        default='unmatch2match_train')
    parser.add_argument('--val_name', type=str, default='unmatch2match_val')
    parser.add_argument('--test_name', type=str, default='unmatch2match_test')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--hyper_search', action='store_false', default=True)
    parser.add_argument('--verbose', action='store_false', default=True)
    parser.add_argument('--covariates',
                        type=List,
                        default=['AGE', 'SEX', 'MMSE', 'DX'])

    args, _ = parser.parse_known_args()
    return args


def covariatess_effects_estimator(args):
    """
    Main function for estimating covariates effects
    """
    ROIs = txt2list(global_config.ROI_features_path)
    cols = txt2list(global_config.columns_path)
    for fold in range(args.nb_folds):
        # run for each fold
        if args.verbose:
            print('Fitting covariates effects estimator for fold', fold)
        fold_input_path = os.path.join(args.data_path, args.dataset_pair,
                                       str(fold))
        fold_output_path = os.path.join(args.output_path, args.model,
                                        args.dataset_pair, str(fold))
        fold_checkpoint_path = os.path.join(args.checkpoint_path, args.model,
                                            args.dataset_pair, str(fold))
        fold_hyper_params_path = os.path.join(args.hyper_params_path,
                                              args.model, args.dataset_pair,
                                              str(fold))
        if args.hyper_search:
            hyper_params = dict()
        else:
            hyper_params = load_pkl(
                os.path.join(fold_hyper_params_path, 'f(x)_hyper_params.pkl'))
        create_folder(fold_output_path)
        create_folder(fold_checkpoint_path)

        # 1. load data
        train_df = pd.read_csv(
            os.path.join(fold_input_path, args.train_name + '.csv'))
        val_df = pd.read_csv(
            os.path.join(fold_input_path, args.val_name + '.csv'))
        test_df = pd.read_csv(
            os.path.join(fold_input_path, args.test_name + '.csv'))
        # 2. data processing
        norm_train_df, norm_val_df, norm_test_df, roi_mean, roi_std = \
            data_proc_estimate_covar_effects(train_df, val_df, test_df,
                                             args.covariates, ROIs)
        # 3. get batch_col for regressing out linear site effect
        batch_col = np.concatenate(
            [norm_train_df['SITE'].values, norm_val_df['SITE'].values], axis=0)
        batch_design = np.vstack([batch_col, np.ones(len(batch_col))]).T
        # 4. placeholder for covariates effects
        train_pred = np.zeros((train_df.shape[0], len(ROIs)))
        val_pred = np.zeros((val_df.shape[0], len(ROIs)))
        test_pred = np.zeros((test_df.shape[0], len(ROIs)))
        # 5. run global checking to get reliable covs
        reliable_covs = global_checker(norm_train_df, norm_val_df, ROIs,
                                       args.covariates, covs_ref_dict)
        # save reliable covariates
        list2txt(reliable_covs,
                 os.path.join(fold_checkpoint_path, 'f(x)_reliable_covs.txt'))
        if len(reliable_covs) == 0:
            if args.verbose:
                print('No reliable covariates!')
            save_covs_effects(fold_output_path, train_df, train_pred, roi_mean,
                              roi_std, ROIs, cols,
                              args.train_name + args.sufix + '.csv')
            save_covs_effects(fold_output_path, val_df, val_pred, roi_mean,
                              roi_std, ROIs, cols,
                              args.val_name + args.sufix + '.csv')
            save_covs_effects(fold_output_path, test_df, test_pred, roi_mean,
                              roi_std, ROIs, cols,
                              args.test_name + args.sufix + '.csv')
        else:
            if args.verbose:
                print('Estimating covariates effects with:', reliable_covs)
            for i, roi in enumerate(ROIs):
                if args.verbose:
                    print('...estimating on', i, roi, datetime.now())
                # 1. regress out linear site effect
                train_roi_desite, val_roi_desite, test_roi_desite =\
                    regress_out_site(norm_train_df, norm_val_df, norm_test_df,
                                     roi, batch_design)
                norm_train_df[roi] = train_roi_desite
                norm_val_df[roi] = val_roi_desite
                norm_test_df[roi] = test_roi_desite
                # 2. train XGBoost for covariates effect estimate
                if args.hyper_search:
                    estimator = XGBWrapper(reliable_covs, roi, 'regress')
                    model, best_hyperparams = estimator.tune(
                        norm_train_df, norm_val_df)
                    hyper_params[roi] = best_hyperparams
                else:
                    estimator = XGBWrapper(reliable_covs,
                                           roi,
                                           'regress',
                                           hyper_params=hyper_params[roi])
                    model = estimator.train(norm_train_df)
                # 3. make prediction
                train_roi_pred = estimator.predict(model, norm_train_df)
                val_roi_pred = estimator.predict(model, norm_val_df)
                test_roi_pred = estimator.predict(model, norm_test_df)
                # 4. reliablility voting
                vote = reliability_voting(val_roi_pred.reshape((-1, )),
                                          val_roi_desite.reshape((-1, )),
                                          corr_type='pearsonr')
                if vote == 1:
                    train_pred[:, i] = train_roi_pred
                    val_pred[:, i] = val_roi_pred
                    test_pred[:, i] = test_roi_pred
                else:
                    # downgrade to linear estimator
                    train_roi_lin, val_roi_lin, test_roi_lin =\
                        linear_covariates_effects_estimator(
                            train_roi_desite, val_roi_desite, test_roi_desite,
                            norm_train_df[reliable_covs].values,
                            norm_val_df[reliable_covs].values,
                            norm_test_df[reliable_covs].values)
                    train_pred[:, i] = train_roi_lin
                    val_pred[:, i] = val_roi_lin
                    test_pred[:, i] = test_roi_lin
            # save prediction
            save_covs_effects(fold_output_path, train_df, train_pred, roi_mean,
                              roi_std, ROIs, cols,
                              args.train_name + args.sufix + '.csv')
            save_covs_effects(fold_output_path, val_df, val_pred, roi_mean,
                              roi_std, ROIs, cols,
                              args.val_name + args.sufix + '.csv')
            save_covs_effects(fold_output_path, test_df, test_pred, roi_mean,
                              roi_std, ROIs, cols,
                              args.test_name + args.sufix + '.csv')
            # save searched hyper_params
            if args.hyper_search:
                save_pkl(
                    hyper_params,
                    os.path.join(fold_checkpoint_path,
                                 'f(x)_hyper_params.pkl'))


if __name__ == '__main__':
    covariatess_effects_estimator(covariates_effects_estimator_parser())
