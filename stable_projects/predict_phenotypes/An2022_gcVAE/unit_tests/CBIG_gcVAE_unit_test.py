#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import torch
import warnings
import unittest
import numpy as np
from model.goalDNN import goalDNN
from model.VAE import VAE
from utils.misc import one_hot
from utils.check_reference_results import\
    check_reference_results, check_results_args_parser
from examples.CBIG_gcVAE_example import example_wrapper, example_args_parser


class CBIG_gcVAE_unit_test(unittest.TestCase):
    def test_goaldnn_init(self):
        """
        Test initialization of goalDNN
        """
        torch.manual_seed(0)
        model = goalDNN(
            in_dim=20,
            nb_category=3,
            nb_measures=1,
            p_dropout=0.1,
            hidden_dims=[64, 32])
        rng = np.random.default_rng(seed=0)
        x = torch.tensor(rng.random(size=(10, 20))).float()
        [mmse_hat, dx_hat] = model(x)
        mmse_hat_mean = torch.mean(mmse_hat).detach().numpy()
        dx_hat_mean = torch.mean(dx_hat).detach().numpy()
        self.assertAlmostEqual(mmse_hat_mean, 0.0588, 4)
        self.assertAlmostEqual(dx_hat_mean, 0.1866, 4)

    def test_cvae_init(self):
        """
        Test initialization of cVAE
        """
        torch.manual_seed(0)
        cvae = VAE(
            in_dim=20,
            nb_classes=2,
            latent_dim=16,
            p_dropout=0.1,
            hidden_dims=[64, 32])
        rng = np.random.default_rng(seed=0)
        x = torch.tensor(rng.random(size=(10, 20))).float()
        rng = np.random.default_rng(seed=0)
        y = rng.integers(low=0, high=2, size=10)
        y_onehot = one_hot(y, nb_classes=2)
        y_onehot = torch.tensor(y_onehot).float()
        [x_hat, _, _, _] = cvae(x, y_onehot, 0)
        x_hat_mean = torch.mean(x_hat).detach().numpy()
        self.assertAlmostEqual(x_hat_mean, -0.0245, 4)

    def test_training(self):
        """
        Test the training of goalDNN, cVAE, gcVAE, XGBoost
        """
        warnings.filterwarnings('ignore')
        example_args = example_args_parser()
        example_args.unittest = True
        # prepare data
        example_args.stage = 'prepare'
        example_wrapper(example_args)
        # train goalDNN model
        example_args.stage = 'train'
        example_args.model = 'goalDNN'
        example_wrapper(example_args)
        # train cVAE model
        example_args.model = 'cVAE'
        example_wrapper(example_args)
        # train gcVAE model
        example_args.model = 'gcVAE'
        example_wrapper(example_args)
        # cVAE harmonization
        example_args.stage = 'predict'
        example_args.model = 'cVAE'
        example_wrapper(example_args)
        # gcVAE harmonization
        example_args.model = 'gcVAE'
        example_wrapper(example_args)
        # goalDNN prediction
        example_args.model = 'goalDNN'
        example_wrapper(example_args)
        # XGBoost
        example_args.stage = 'train'
        example_args.model = 'XGBoost'
        example_wrapper(example_args)
        # compare with reference results
        check_args = check_results_args_parser()
        check_args.unittest = True
        check_reference_results(check_args)


if __name__ == '__main__':
    unittest.main()
