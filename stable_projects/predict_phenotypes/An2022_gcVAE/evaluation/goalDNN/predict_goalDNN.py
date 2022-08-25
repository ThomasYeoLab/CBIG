#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import argparse
import torch
from config import global_config
from utils.misc import load_pkl, txt2list, save_pkl, create_folder
from utils.nn_misc import extract_goaldnn_input


def predict_goalDNN_args_parser():
    parser = argparse.ArgumentParser(prog='PredgoalDNNArgs')
    # parameteres
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--checkpoint_path', type=str, default='/')
    parser.add_argument('--save_path', type=str, default='/')
    parser.add_argument('--exp', type=str, default='sample_size')
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--fold', type=int, default=0)
    parser.add_argument(
        '--harm_models', type=list, default=['ComBat', 'cVAE', 'gcVAE'])
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')
    parser.add_argument('--testsets', type=list, default=[])

    args, _ = parser.parse_known_args()
    return args


def predict(args):
    """
    Making prediction using trained goalDNN model

    Args:
        args (tuple): Parameters
    """
    # load model
    model = torch.load(
        os.path.join(args.checkpoint_path, str(args.fold), 'goalDNN.pt'),
        map_location='cpu')
    model.to(torch.device('cpu'))
    model.eval()
    # read features of training data ==> 108 brain ROI volumes
    ROI_features = txt2list(global_config.ROI_features_path)
    for testset in args.testsets:
        # we have multiple testsets to make prediction
        test_pkl = load_pkl(
            os.path.join(args.data_path, str(args.fold),
                         'test_' + testset + '.pkl'))
        mean = test_pkl['mean']
        std = test_pkl['std']
        testROIs, testMMSEs, testDXs, _ = \
            extract_goaldnn_input(test_pkl, ROI_features)
        # making prediction
        [testMMSEs_pred, testDXs_pred] = model(testROIs.float())
        testMMSEs = testMMSEs.numpy()
        testMMSEs_pred = testMMSEs_pred.detach().numpy()
        testDXs = testDXs.numpy()
        testDXs_pred = testDXs_pred.detach().numpy()
        # save prediction
        test_pred = dict()
        test_pred['MMSE'] = dict()
        test_pred['MMSE']['Pred'] = testMMSEs_pred * std['MMSE'] + mean['MMSE']
        test_pred['MMSE']['Pred'][test_pred['MMSE']['Pred'] < 0] = 0
        test_pred['MMSE']['Pred'][test_pred['MMSE']['Pred'] > 30] = 30
        test_pred['MMSE']['GT'] = testMMSEs * std['MMSE'] + mean['MMSE']
        test_pred['DX'] = dict()
        test_pred['DX']['Pred'] = testDXs_pred
        test_pred['DX']['GT'] = testDXs
        test_pred['RID'] = test_pkl['RID']
        create_folder(os.path.join(args.save_path, str(args.fold)))
        # save to pkl file
        save_pkl(
            test_pred,
            os.path.join(args.save_path, str(args.fold),
                         'pred_test_' + testset + '.pkl'))


def predict_wrapper(args):
    """
    Wrapper function for making prediction

    Args:
        args (tuple): Parameters
    """
    if args.exp == 'sample_size':
        args.harm_models = ['ComBat', 'cVAE', 'gcVAE']
    else:
        args.harm_models = ['ComBat', 'ComBat4cov', 'cVAE', 'gcVAE']
    nonADNI_dataset = args.dataset_pair.split('-')[1]
    testsets = ['ADNI_unharm', nonADNI_dataset + '_unharm']
    for harm_model in args.harm_models:
        testsets.append(nonADNI_dataset + harm_model + '_harm')
    args.testsets = testsets
    for fold in range(args.nb_folds):
        args.fold = fold
        predict(args)


if __name__ == '__main__':
    predict_wrapper(predict_goalDNN_args_parser())
