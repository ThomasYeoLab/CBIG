#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import torch
import warnings
import unittest
import numpy as np
import pandas as pd
from model.VAE import VAE
from utils.misc import one_hot
from examples.DeepResBat_example import example_args_parser, example_wrapper


class CBIG_DeepResBat_unit_test(unittest.TestCase):

    def test_cvae_init(self):
        """
        Test initialization of cVAE
        """
        torch.manual_seed(0)
        cvae = VAE(in_dim=20,
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

    def test_deepresbat(self):
        """
        Test the DeepResBat model
        """
        warnings.filterwarnings('ignore')
        example_args = example_args_parser()
        example_args.unittest = True
        example_args.gen_data = True
        example_args.covar_effects = True
        example_args.gen_res = True
        example_args.harm_res = True
        example_args.add_back = True
        example_args.manova = True
        example_wrapper(example_args)
        ref_results = np.array([0.416988, 0.455762, 0.625207, 0.499431])
        manova_results_path = os.path.join(example_args.working_dir, 'results',
                                           'assoc_manova')
        manova_csv_path = os.path.join(manova_results_path, 'MMSE_demean',
                                       example_args.dataset_pair,
                                       'p_DeepResBat_ADNI-AIBL.csv')
        manova_df = pd.read_csv(manova_csv_path)
        results = manova_df['PillarTrace'].values.reshape((-1, ))
        self.assertTrue(np.allclose(ref_results, results, atol=0.1))


if __name__ == '__main__':
    unittest.main()
