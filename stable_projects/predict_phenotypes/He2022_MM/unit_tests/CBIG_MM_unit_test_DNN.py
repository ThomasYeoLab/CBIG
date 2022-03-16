#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import glob
import shutil
import unittest
import subprocess
import numpy as np

if 'CBIG_CODE_DIR' in os.environ:
    ut_dir = os.path.join(os.environ['CBIG_CODE_DIR'], 'stable_projects',
                          'predict_phenotypes', 'He2022_MM', 'unit_tests')
else:
    ut_dir = os.getcwd()


class TestMM(unittest.TestCase):
    def test_DNN_MM_training(self):
        input_dir = os.path.join(ut_dir, 'data')
        output_dir = os.path.join(ut_dir, 'output', 'output_dnn')
        subprocess.call([
            'python3', '../cbig/He2022/CBIG_ukbb_dnn.py', '--exp-dataset',
            '--in_dir', input_dir, '--out_dir', output_dir, '--inter_dir',
            input_dir, '--gpu', '0', '--seed', '1', '--epochs', '30',
            '--metric', 'cod', '--weight_decay', '8.447320e-04', '--lr',
            '3.645653e-03', '--dropout', '0.241964', '--scheduler_decrease',
            '312', '--n_l1', '87', '--n_l2', '386', '--n_l3', '313',
            '--n_layer', '3'
        ])
        ref = os.path.join(ut_dir, 'results', 'dnn_exp_dataset_base.npz')
        res = os.path.join(output_dir, 'dnn_exp_dataset_base.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(res)[0])
        ref_val = np.mean(ref['val_cor_record'])
        res_val = np.mean(res['val_cor_record'])
        self.assertTrue(
            np.isclose(ref_val, res_val),
            'DNN base training metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'output'))

    def test_DNN_MM(self):
        input_dir = os.path.join(ut_dir, 'data')
        output_dir = os.path.join(ut_dir, 'output', 'output_dnn')
        os.makedirs(output_dir, exist_ok=True)
        shutil.copy(
            os.path.join(input_dir, 'dnn_exp_dataset_base.npz'),
            os.path.join(output_dir, 'dnn_exp_dataset_base.npz'))
        subprocess.call([
            'python3', '../cbig/He2022/CBIG_ukbb_dnn_mm.py', '--exp-dataset',
            '--in_dir', input_dir, '--out_dir', output_dir, '--inter_dir',
            input_dir, '--rng', '2'
        ])
        ref = os.path.join(ut_dir, 'results',
                           'meta_result_test_exp_dataset.npz')
        res = os.path.join(output_dir, 'meta_result_test_exp_dataset.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(res)[0])
        ref_val = np.mean(ref['meta_cor'])
        res_val = np.mean(res['meta_cor'])
        self.assertTrue(
            np.isclose(ref_val, res_val),
            'DNN MM metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'output'))

    def test_DNN_stacking(self):
        input_dir = os.path.join(ut_dir, 'data')
        output_dir = os.path.join(ut_dir, 'output', 'output_dnn')
        os.makedirs(output_dir, exist_ok=True)
        shutil.copy(
            os.path.join(input_dir, 'dnn_exp_dataset_base.npz'),
            os.path.join(output_dir, 'dnn_exp_dataset_base.npz'))
        subprocess.call([
            'python3', '../cbig/He2022/CBIG_ukbb_dnn_mm_stacking.py',
            '--exp-dataset', '--in_dir', input_dir, '--out_dir', output_dir,
            '--inter_dir', input_dir, '--rng', '2'
        ])
        ref = os.path.join(ut_dir, 'results', 'meta_stacking_result_test.npz')
        res = os.path.join(output_dir, 'meta_stacking_result_test.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(res)[0])
        ref_val = np.mean(ref['meta_cor'])
        res_val = np.mean(res['meta_cor'])
        self.assertTrue(
            np.isclose(ref_val, res_val),
            'DNN MM stacking metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'output'))


if __name__ == '__main__':
    unittest.main()
