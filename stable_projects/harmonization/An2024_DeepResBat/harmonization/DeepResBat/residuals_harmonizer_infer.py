#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import argparse
from harmonization.cVAE.predict_cVAE import predict_cvae_args_parser, predict


def residuals_harmonizer_infer_args_parser():
    """
    Parameters for making prediction using trained cVAE model
    """
    parser = argparse.ArgumentParser(prog='PredcVAEArgs')
    # input and output path
    parser.add_argument('--raw_data_path', type=str, default='/')
    parser.add_argument('--harm_input_path', type=str, default='/')
    parser.add_argument('--checkpoint_path', type=str, default='/')
    parser.add_argument('--harm_output_path', type=str, default='/')
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')
    parser.add_argument('--exp', type=str, default='unmatch2match')
    parser.add_argument('--nb_pred', type=int, default=100)
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--model_name', type=str, default='DeepResBat')
    parser.add_argument('--harmonizer', type=str, default='cVAE')
    parser.add_argument('--sufix', type=str, default='_G')
    parser.add_argument('--harm_files',
                        type=list,
                        default=['train', 'val', 'test'])
    parser.add_argument('--origin_files',
                        type=list,
                        default=['train', 'val', 'test'])
    # in case there are unexcepted args
    pred_args, _ = parser.parse_known_args()

    return pred_args


def residuals_harmonizer_infer(args):
    """
    Using trained residual harmoinzer to infer
    By deafult we choose cVAE as residual harmonization model
    """
    if args.harmonizer == 'cVAE':
        train_cvae_args = predict_cvae_args_parser()
        # fit all parameters
        for arg in vars(train_cvae_args):
            setattr(train_cvae_args, arg, getattr(args, arg))
        # train cVAE model
        predict(train_cvae_args)
    else:
        raise NotImplementedError(
            "Only support cVAE as residual harmonizer for now")


if __name__ == '__main__':
    residuals_harmonizer_infer(residuals_harmonizer_infer_args_parser())
