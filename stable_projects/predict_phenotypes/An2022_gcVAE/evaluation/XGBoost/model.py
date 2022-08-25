#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import xgboost as xgb
from utils.metrics import site_prediction_metric


def site_pred_model(args, train_df, val_df, ROIs):
    """
    XGBoost model for predicting site

    Args:
        args (tuple): Parameters
        train_df (class DataFrame): Dataframe for training model
        val_df (class DataFrame): Dataframe for validation
        ROIs (list): Features
    """
    # training set
    x_train = train_df[ROIs]
    y_train = train_df['SITE']
    xg_train = xgb.DMatrix(x_train.values, y_train.values, feature_names=ROIs)
    # validation test
    x_val = val_df[ROIs]
    y_val = val_df['SITE']
    xg_val = xgb.DMatrix(x_val.values, y_val.values, feature_names=ROIs)

    # train XGBoost
    model = xgb.train(
        params=args.params,
        dtrain=xg_train,
        num_boost_round=args.num_boost_rounds,
        verbose_eval=False)

    # making prediction
    val_pred = model.predict(xg_val)
    # evaluation
    val_truth = y_val.values
    opt_score = 0
    opt_threshold = 0
    for threshold in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5):
        _, score = site_prediction_metric(val_df['RID'].values, val_pred,
                                          val_truth, threshold)
        if score > opt_score:
            opt_score = score
            opt_threshold = threshold

    return opt_score, opt_threshold, model
