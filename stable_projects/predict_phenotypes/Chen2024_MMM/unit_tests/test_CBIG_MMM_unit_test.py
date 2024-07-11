#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import sys
import shutil
import unittest
import subprocess
import numpy as np
sys.path.append('../cbig/Chen2024')
from config import config
from CBIG_misc import compare_dicts, get_phe_num

ut_dir = config.UT_DIR
in_dir = config.IN_DIR_UT
out_dir = config.OUT_DIR_UT
inter_dir = config.INTER_DIR_UT
model_dir = config.MODEL_DIR_UT
ref_dir = config.REF_OUT_DIR_UT

class TestMMM(unittest.TestCase):
    def test_DNN(self):
        subprocess.call([
            'python', '../cbig/Chen2024/CBIG_dnn_xlarge_train.py', '--exp-dataset',
            '--unit-test', '--seed', '1', '--epochs', '100', '--metric', 'cod',
            '--weight_decay', '1e-7', '--lr', '0.01', '--dropout', '0.2',
            '--n_l1', '64', '--n_l2', '64', '--n_l3', '64', '--n_hidden_layer', '2',
            '--batch_size', '64', '--patience', '25'
        ])
        ref = np.load(os.path.join(ref_dir, 'output_intermediate', 'dnn_base.npz'), allow_pickle=True)
        res = np.load(os.path.join(inter_dir, 'dnn_base.npz'), allow_pickle=True)
        self.assertTrue(
            compare_dicts(ref, res),
            'DNN base training metric ' + str(res) + ' vs. ref ' + str(ref))

        subprocess.call([
            'python', '../cbig/Chen2024/CBIG_dnn_xlarge_predict.py',
            '--exp-dataset', '--unit-test'
        ])
        ref = np.load(os.path.join(ref_dir, 'output_intermediate', 'dnn_prediction.npz'))
        res = np.load(os.path.join(inter_dir, 'dnn_prediction.npz'))
        self.assertTrue(
            compare_dicts(ref, res),
            'DNN prediction ' + str(res) + ' vs. ref ' + str(ref))

        shutil.copyfile(os.path.join(in_dir, 'exp_test_split_ind_3.npz'),
                        os.path.join(inter_dir, 'exp_test_split_ind_3.npz'))
        subprocess.call([
            'python', '../cbig/Chen2024/CBIG_dnn_transfer_learning.py',
            '--exp-dataset', '--unit-test'
        ])
        ref = np.load(os.path.join(ref_dir, 'transfer_learning_2exp_test_result.npz'), allow_pickle=True)
        res = np.load(os.path.join(out_dir, 'transfer_learning_2exp_test_result.npz'), allow_pickle=True)
        self.assertTrue(
            compare_dicts(ref, res),
            'Transfer learning results ' + str(res) + ' vs. ref ' + str(ref))

    def test_stacking(self):
        os.makedirs(inter_dir, exist_ok=True)
        shutil.copyfile(os.path.join(ref_dir, 'output_intermediate', 'dnn_prediction.npz'),
                        os.path.join(inter_dir, 'dnn_prediction.npz'))
        shutil.copyfile(os.path.join(ref_dir, 'output_intermediate', 'rr_prediction.npz'),
                        os.path.join(inter_dir, 'rr_prediction.npz'))
        shutil.copyfile(os.path.join(in_dir, 'exp_test_split_ind_3.npz'),
                        os.path.join(inter_dir, 'exp_test_split_ind_3.npz'))

        n_phe_L = get_phe_num(in_dir, 'exp_train_L')
        n_phe_M = get_phe_num(in_dir, 'exp_train_M')

        for phe_idx in range(n_phe_L):
            subprocess.call([
                'python', '../cbig/Chen2024/CBIG_rr_large.py', '--exp-dataset',
                '--unit-test', '--phe_idx', str(phe_idx)
            ])
        for phe_idx in range(n_phe_M):
            subprocess.call([
                'python', '../cbig/Chen2024/CBIG_rr_medium.py', '--exp-dataset',
                '--unit-test', '--phe_idx', str(phe_idx)
            ])
        for phe_idx in range(n_phe_L):
            subprocess.call([
                'python', '../cbig/Chen2024/CBIG_rr_large.py', '--exp-dataset',
                '--unit-test', '--phe_idx', str(phe_idx), '--2layer'
            ])
        for phe_idx in range(n_phe_M):
            subprocess.call([
                'python', '../cbig/Chen2024/CBIG_rr_medium.py', '--exp-dataset',
                '--unit-test', '--phe_idx', str(phe_idx), '--2layer'
            ])

        subprocess.call([
            'python', '../cbig/Chen2024/CBIG_mm_stacking.py', '--exp-dataset',
            '--unit-test', '--log_stem', 'MM_stacking'
        ])

        ref = np.load(os.path.join(ref_dir, 'MM_stacking_2exp_test_result.npz'), allow_pickle=True)
        res = np.load(os.path.join(out_dir, 'MM_stacking_2exp_test_result.npz'), allow_pickle=True)
        self.assertTrue(
            compare_dicts(ref, res),
            'Meta-matching with stacking results ' + str(res) + ' vs. ref ' + str(ref))

        subprocess.call([
            'python', '../cbig/Chen2024/CBIG_mm_stacking.py', '--exp-dataset',
            '--unit-test', '--log_stem', 'dataset_stacking'
        ])

        ref = np.load(os.path.join(ref_dir, 'dataset_stacking_2exp_test_result.npz'), allow_pickle=True)
        res = np.load(os.path.join(out_dir, 'dataset_stacking_2exp_test_result.npz'), allow_pickle=True)
        self.assertTrue(
            compare_dicts(ref, res),
            'Meta-matching with dataset stacking results ' + str(res) + ' vs. ref ' + str(ref))

        subprocess.call([
            'python', '../cbig/Chen2024/CBIG_mm_stacking.py', '--exp-dataset',
            '--unit-test', '--log_stem', 'multilayer_stacking'
        ])

        ref = np.load(os.path.join(ref_dir, 'multilayer_stacking_2exp_test_result.npz'), allow_pickle=True)
        res = np.load(os.path.join(out_dir, 'multilayer_stacking_2exp_test_result.npz'), allow_pickle=True)
        self.assertTrue(
            compare_dicts(ref, res),
            'Multilayer meta-matching results ' + str(res) + ' vs. ref ' + str(ref))

    def tearDown(self):
        shutil.rmtree(model_dir)
        shutil.rmtree(inter_dir)
        shutil.rmtree(out_dir)

if __name__ == '__main__':
    unittest.main()