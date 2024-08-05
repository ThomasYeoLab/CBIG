#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import argparse
from harmonization.cVAE.train_cVAE import train_cvae_args_parser, train


def residuals_harmonizer_fit_args_parser():
    """
    Parameters for training harmonizer (cVAE model)
    """
    parser = argparse.ArgumentParser(prog='TrueBatGArgs')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--GPU', type=int, default=1)
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--checkpoint_path', type=str, default='/')
    parser.add_argument('--isSaving', action='store_true', default=False)
    parser.add_argument('--model_name', type=str, default='DeepResBat')
    parser.add_argument('--harmonizer', type=str, default='cVAE')
    parser.add_argument('--sufix', type=str, default='_G')
    parser.add_argument('--cpu', action='store_true', default=False)
    # training parameters
    parser.add_argument('--epochs', type=int, default=1000)
    parser.add_argument('--adv_epochs', type=int, default=1)
    parser.add_argument('--step', type=int, default=0)
    parser.add_argument('--batch_size', type=int, default=512)
    parser.add_argument('--nb_classes', type=int, default=2)
    parser.add_argument('--in_dim', type=int, default=108)
    parser.add_argument('--VAE_hidden_dims', type=list, default=[96, 64])
    parser.add_argument('--D_hidden_dims', type=list, default=[32, 32])
    parser.add_argument('--dis_lr', type=float, default=1e-2)
    parser.add_argument('--coff', type=float, default=0)
    # hyper-parameters
    parser.add_argument('--lr', type=float, default=1e-2)
    parser.add_argument('--drop_out', type=float, default=0.1)
    parser.add_argument('--alpha', type=float, default=0.001)
    parser.add_argument('--lambda_', type=float, default=0.001)
    parser.add_argument('--gamma', type=float, default=1.)
    parser.add_argument('--lr_step', type=int, default=100)
    parser.add_argument('--latent_dim', type=int, default=32)
    parser.add_argument('--h1', type=int, default=512)
    parser.add_argument('--h2', type=int, default=512)
    parser.add_argument('--h3', type=int, default=512)
    parser.add_argument('--h4', type=int, default=512)
    parser.add_argument('--nb_layers', type=int, default=3)

    train_args, _ = parser.parse_known_args()
    return train_args


def residuals_harmonizer_fit(args):
    """
    Fitting residual harmoinzer
    By deafult we chose cVAE as residual harmonization model
    """
    if args.harmonizer == 'cVAE':
        train_cvae_args = train_cvae_args_parser()
        # fit all parameters
        for arg in vars(train_cvae_args):
            setattr(train_cvae_args, arg, getattr(args, arg))
        # train cVAE model
        if args.isSaving:
            train(train_cvae_args)
        else:
            return train(train_cvae_args)
    else:
        raise NotImplementedError(
            "Only support cVAE as residual harmonizer for now")


if __name__ == '__main__':
    residuals_harmonizer_fit(residuals_harmonizer_fit_args_parser())
