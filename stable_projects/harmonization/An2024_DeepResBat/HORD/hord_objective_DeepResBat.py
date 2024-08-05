#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import torch
import argparse
import pandas as pd
import numpy as np
from config import global_config
from harmonization.DeepResBat.residuals_harmonizer_fit import\
    residuals_harmonizer_fit, residuals_harmonizer_fit_args_parser
from model.XGB import XGBWrapper
from utils.normalization import ICV_norm
from utils.misc import load_pkl, txt2list, one_hot
from utils.nn_misc import vae_harm_predict
from utils.metrics import data_pred_metric


def objHORD_args_parser():
    """
    [Feeding parameters]
    """
    # general parameters
    objHORD_parser = argparse.ArgumentParser(prog='objHORDArgs')
    objHORD_parser.add_argument('--data_path', type=str, required=True)
    objHORD_parser.add_argument('--harm_output_path', type=str, required=True)
    objHORD_parser.add_argument('--ysf_sufix', type=str, required=True)
    objHORD_parser.add_argument('--GPU', type=int, required=True)
    # hyperparameters for cVAE
    objHORD_parser.add_argument('--lr', type=float, required=True)
    objHORD_parser.add_argument('--drop_out', type=float, required=True)
    objHORD_parser.add_argument('--alpha', type=float, required=True)
    objHORD_parser.add_argument('--lambda_', type=float, required=True)
    objHORD_parser.add_argument('--gamma', type=float, required=True)
    objHORD_parser.add_argument('--lr_step', type=int, required=True)
    objHORD_parser.add_argument('--latent_dim', type=int, required=True)
    objHORD_parser.add_argument('--h1', type=int, required=True)
    objHORD_parser.add_argument('--h2', type=int, required=True)
    objHORD_parser.add_argument('--h3', type=int, required=True)
    objHORD_parser.add_argument('--h4', type=int, required=True)
    objHORD_parser.add_argument('--nb_layers', type=int, default=3)

    objHORD_args, _ = objHORD_parser.parse_known_args()
    return objHORD_args


def objHORD(objHORD_args):
    """
    [cVAE objective function for HORD optimization]

    Args:
        objHORD_args ([type]): [description]
    """
    # step0, data preparation
    ROIs = txt2list(global_config.ROI_features_path)
    train_pkl = load_pkl(
        os.path.join(objHORD_args.data_path,
                     'unmatch2match_train' + objHORD_args.ysf_sufix + '.pkl'))
    val_pkl = load_pkl(
        os.path.join(objHORD_args.data_path,
                     'unmatch2match_val' + objHORD_args.ysf_sufix + '.pkl'))
    mean = train_pkl['mean'][ROIs].values
    std = train_pkl['std'][ROIs].values
    # train
    train_array = train_pkl['data'][ROIs].values
    train_site_onehot = one_hot(train_pkl['data']['SITE'].values)
    train_tensor = np.concatenate((train_array, train_site_onehot), axis=1)
    train_tensor = torch.tensor(train_tensor).float()
    # val
    val_array = val_pkl['data'][ROIs].values
    val_site_onehot = one_hot(val_pkl['data']['SITE'].values)
    val_tensor = np.concatenate((val_array, val_site_onehot), axis=1)
    val_tensor = torch.tensor(val_tensor).float()
    # step1, train cVAE
    train_cVAE_args = residuals_harmonizer_fit_args_parser()
    train_cVAE_args.data_path = objHORD_args.data_path
    if global_config.gpuQ:
        train_cVAE_args.GPU = -1
    else:
        train_cVAE_args.GPU = objHORD_args.GPU
    # parse hyperparameters
    train_cVAE_args.sufix = objHORD_args.ysf_sufix
    train_cVAE_args.lr = objHORD_args.lr
    train_cVAE_args.drop_out = objHORD_args.drop_out
    train_cVAE_args.alpha = objHORD_args.alpha
    train_cVAE_args.lambda_ = objHORD_args.lambda_
    train_cVAE_args.gamma = objHORD_args.gamma
    train_cVAE_args.lr_step = objHORD_args.lr_step
    train_cVAE_args.latent_dim = objHORD_args.latent_dim
    train_cVAE_args.h1 = objHORD_args.h1
    train_cVAE_args.h2 = objHORD_args.h2
    train_cVAE_args.h3 = objHORD_args.h3
    train_cVAE_args.h4 = objHORD_args.h4
    train_cVAE_args.nb_layers = objHORD_args.nb_layers
    # get trained cVAE model & valROIMSE
    cVAE, valROIMSE = residuals_harmonizer_fit(train_cVAE_args)
    device = torch.device('cpu')
    cVAE = cVAE.to(device)  # move to CPU
    # step2, making prediction
    train_pred = np.zeros_like(train_array)
    val_pred = np.zeros_like(val_array)
    for i in range(10):
        train_pred += vae_harm_predict(cVAE, train_tensor, i) / 10.0
        val_pred += vae_harm_predict(cVAE, val_tensor, i) / 10.0
    train_pred = (train_pred * std) + mean
    val_pred = (val_pred * std) + mean
    # step2, data preparation for XGBoost site prediction
    train_df_g = pd.DataFrame(data=train_pred, columns=ROIs)
    train_df_g['RID'] = train_pkl['RID'].values
    train_df_g['SITE'] = train_pkl['data']['SITE'].values
    val_df_g = pd.DataFrame(data=val_pred, columns=ROIs)
    val_df_g['SITE'] = val_pkl['data']['SITE'].values
    val_df_g['RID'] = val_pkl['RID'].values
    # run ICV normalization
    norm_train_df = ICV_norm(train_df_g, ROIs)
    norm_val_df = ICV_norm(val_df_g, ROIs)
    # step3, run XGBoost site prediction
    predictor = XGBWrapper(features=ROIs,
                           target='SITE',
                           task='site',
                           nb_boost=100)
    model, best_hyper_params = predictor.tune(norm_train_df, norm_val_df)
    val_pred = predictor.predict(model, norm_val_df)
    _, valSiteAcc = data_pred_metric(norm_val_df['RID'].values, val_pred,
                                     norm_val_df['SITE'].values,
                                     best_hyper_params['threshold'])
    # output result
    print(valROIMSE, valSiteAcc, sep=',')


if __name__ == '__main__':
    objHORD(objHORD_args_parser())
