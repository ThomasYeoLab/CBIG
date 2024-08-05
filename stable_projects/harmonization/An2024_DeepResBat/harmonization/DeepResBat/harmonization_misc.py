#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import copy
import random
import numpy as np
from scipy import stats
from model.XGB import XGBWrapper
from utils.misc import save_df
from utils.imputation import ff_multi_cols
from utils.normalization import z_norm


def reliability_voting(pred,
                       target,
                       threshold=0.95,
                       ratio=0.5,
                       corr_type='pearsonr',
                       nb_votes=100):
    """
    A voting function to check whether prediction is reliable or not
    """
    votes = []
    nb_samples = pred.shape[0]
    for i in range(nb_votes):
        indexs = [i for i in range(nb_samples)]
        random.Random(i).shuffle(indexs)
        selected_indexs = indexs[:int(nb_samples * ratio)]
        if corr_type == 'pearsonr':
            _, p = stats.pearsonr(pred[selected_indexs],
                                  target[selected_indexs])
        elif corr_type == 'spearmanr':
            _, p = stats.spearmanr(pred[selected_indexs],
                                   target[selected_indexs])
        else:
            raise NotImplementedError("Only support pearsonr & spearmanr!")
        if p >= 0.05:
            votes.append(0)
        else:
            votes.append(1)
    if np.sum(votes) / len(votes) >= threshold:
        return 1
    else:
        return 0


def global_checker(norm_train_df, norm_val_df, ROIs, covariates, cov_ref_dict):
    """
    Use all ROIs to predict covariates to get globally reliable covariates
    """
    cov_type_dict = {
        "SEX": 'spearmanr',
        "DX": 'spearmanr',
        "AGE": 'pearsonr',
        "MMSE": 'pearsonr'
    }
    reliable_covs = []
    for cov in covariates:
        predictor = XGBWrapper(ROIs, cov, cov_ref_dict[cov])
        # grid search on validation set
        best_model, best_hyperparams = predictor.tune(norm_train_df,
                                                      norm_val_df)
        val_pred = predictor.predict(best_model, norm_val_df)
        if cov_ref_dict[cov] == 'binary':
            best_threshold = best_hyperparams['threshold']
            val_pred = (val_pred >= best_threshold) * 1
        elif cov_ref_dict[cov] == 'multi':
            val_pred = val_pred.argmax(axis=1)
        else:
            pass
        pred_vec = val_pred.reshape((-1, ))
        target_vec = (norm_val_df[cov].values).reshape((-1, ))
        vote = reliability_voting(pred_vec,
                                  target_vec,
                                  corr_type=cov_type_dict[cov])
        if vote == 1:
            reliable_covs.append(cov)
    return reliable_covs


def save_covs_effects(output_path, df, roi_predicted, roi_mean, roi_std, ROIs,
                      cols, save_name):
    """
    Save the covairiates effects
    """
    df = df[cols]
    df[ROIs] = roi_predicted
    df[ROIs] = (df[ROIs] * roi_std) + roi_mean
    save_df(df, os.path.join(output_path, save_name))


def regress_out_site(train_df, val_df, test_df, roi, batch_design):
    """
    Regress out linear site effects
    """
    train_roi_vec = (train_df[roi].values).reshape((-1, 1))
    val_roi_vec = (val_df[roi].values).reshape((-1, 1))
    test_roi_vec = (test_df[roi].values).reshape((-1, 1))
    # estimate beta
    roi_vec = np.concatenate([train_roi_vec, val_roi_vec], axis=0)
    beta, _ = np.linalg.lstsq(batch_design, roi_vec, rcond=None)[0]
    # desite
    train_site, val_site, test_site =\
        train_df['SITE'].values, val_df['SITE'].values, test_df['SITE'].values
    train_roi_desite = train_roi_vec - beta * train_site.reshape((-1, 1))
    val_roi_desite = val_roi_vec - beta * val_site.reshape((-1, 1))
    test_roi_desite = test_roi_vec - beta * test_site.reshape((-1, 1))

    return train_roi_desite, val_roi_desite, test_roi_desite


def linear_covariates_effects_estimator(train_roi, val_roi, test_roi,
                                        train_covs, val_covs, test_covs):
    """
    Incase XGBoost estimator is not reliable,
    Fit a linear model to estimate effects of covariates
    """
    # reshape
    train_roi, val_roi, test_roi = train_roi.reshape((-1, 1)), val_roi.reshape(
        (-1, 1)), test_roi.reshape((-1, 1))
    # fit model with train and validation data
    train_val_roi = np.concatenate([train_roi, val_roi], axis=0)
    train_val_covs = np.concatenate([train_covs, val_covs], axis=0)
    nb_fit_samples = train_val_roi.shape[0]
    nb_covs = train_covs.shape[1]
    train_val_design = np.concatenate(
        ([train_val_covs, np.ones((nb_fit_samples, 1))]), axis=1)
    # fit a linear model
    betas = np.linalg.lstsq(train_val_design, train_val_roi, rcond=None)[0]
    assert betas.shape[0] == (nb_covs + 1), 'Wrong betas!'
    cov_betas = betas[:nb_covs].reshape((1, -1))
    # get prediction
    train_pred = np.sum(train_covs * cov_betas, axis=1).reshape((-1, ))
    val_pred = np.sum(val_covs * cov_betas, axis=1).reshape((-1, ))
    test_pred = np.sum(test_covs * cov_betas, axis=1).reshape((-1, ))

    return train_pred, val_pred, test_pred


def data_proc_estimate_covar_effects(train_df,
                                     val_df,
                                     test_df,
                                     covariates,
                                     ROIs):
    """
    Data processing to estimate covariates effects
    """
    im_train_df = ff_multi_cols(copy.deepcopy(train_df), covariates)
    im_val_df = ff_multi_cols(copy.deepcopy(val_df), covariates)
    im_test_df = ff_multi_cols(copy.deepcopy(test_df), covariates)
    # set sex as [0, 1]
    im_train_df['SEX'] -= 1
    im_val_df['SEX'] -= 1
    im_test_df['SEX'] -= 1
    # 3. z-normalize ROIs and fillna as 0
    roi_mean = im_train_df[ROIs].mean()
    roi_std = im_train_df[ROIs].std()
    # train
    norm_train_df = copy.deepcopy(im_train_df)
    norm_train_df[ROIs] = z_norm(im_train_df[ROIs], roi_mean, roi_std)
    norm_train_df.fillna(value=0, inplace=True)
    # val
    norm_val_df = copy.deepcopy(im_val_df)
    norm_val_df[ROIs] = z_norm(im_val_df[ROIs], roi_mean, roi_std)
    norm_val_df.fillna(value=0, inplace=True)
    # test
    norm_test_df = copy.deepcopy(im_test_df)
    norm_test_df[ROIs] = z_norm(im_test_df[ROIs], roi_mean, roi_std)
    norm_test_df.fillna(value=0, inplace=True)

    return norm_train_df, norm_val_df, norm_test_df, roi_mean, roi_std
