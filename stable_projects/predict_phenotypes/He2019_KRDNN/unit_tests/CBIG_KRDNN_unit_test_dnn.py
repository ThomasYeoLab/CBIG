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

from cbig.He2019.config import config
from cbig.He2019.CBIG_prepare_data import data_ukbb_fnn, data_ukbb_fnn_sex
from cbig.He2019.CBIG_prepare_data import data_ukbb_brainnetcnn, data_hcp_fnn
from cbig.He2019.CBIG_prepare_data import data_ukbb_gcnn, data_ukbb_gcnn_sex
from cbig.He2019.CBIG_prepare_data import data_hcp_brainnetcnn, data_hcp_gcnn
from cbig.He2019.CBIG_prepare_data import data_ukbb_brainnetcnn_sex
from cbig.He2019.CBIG_prepare_data import get_gcnn_graph

if 'CBIG_CODE_DIR' in os.environ:
    ut_dir = os.path.join(os.environ['CBIG_CODE_DIR'], 'stable_projects',
                          'predict_phenotypes', 'He2019_KRDNN', 'unit_tests')
else:
    ut_dir = os.getcwd()


def rreplace(s, old, new, occurrence=1):
    li = s.rsplit(old, occurrence)
    return new.join(li)


def compare_npz(npz, npz_ref, test):
    test.assertTrue(os.path.isfile(npz_ref), 'no reference file for' + npz)
    np1 = np.load(npz)
    np2 = np.load(npz_ref)
    test.assertTrue(
        np.array_equal(np1.files, np2.files), 'wrong files from' + npz)
    for file in np1.files:
        test.assertTrue(
            np.allclose(np1[file], np2[file]),
            'wrong value from' + npz + 'for file' + file)


def compare_npz_dir(npz_dir, test):
    npzs = glob.glob(os.path.join(npz_dir, '*.npz'))
    for npz in npzs:
        npz_ref = rreplace(npz, '/data/', '/results/')
        compare_npz(npz, npz_ref, test)


def compare_graph(graph, graph_ref, test):
    test.assertTrue(os.path.isfile(graph_ref), 'no reference file for' + graph)
    f1 = open(graph, 'r')
    f2 = open(graph_ref, 'r')
    for l1, l2 in zip(f1, f2):
        n1 = np.fromstring(l1, count=3, sep='\t')
        n2 = np.fromstring(l2, count=3, sep='\t')
        test.assertTrue(
            np.allclose(n1, n2), 'wrong value from' + graph + 'for line' + l1)
    f1.close()
    f2.close()


def compare_graph_dir(graph_dir, test):
    graphs = glob.glob(os.path.join(graph_dir, '*.cites'))
    for graph in graphs:
        graph_ref = rreplace(graph, '/data/', '/results/')
        compare_graph(graph, graph_ref, test)


class TestData(unittest.TestCase):
    def test_data_nilearn(self):
        n_folds = config.EXAMPLE_N_FOLDS
        cur_dir = os.path.join(ut_dir, 'data')
        nilearn_npz = os.path.join(cur_dir, 'nilearn_adhd_fc_beh.npz')

        npz = np.load(nilearn_npz)
        corr_mat = npz['corr_mat']
        beh = npz['beh']
        nilearn_tvt_dir = os.path.join(cur_dir,
                                       'tvt')  # training validation testing
        nilearn_cv_dir = os.path.join(cur_dir, 'cv')  # cross validation
        graph_dir = os.path.join(cur_dir, 'graph')

        os.makedirs(nilearn_tvt_dir, exist_ok=True)
        data_ukbb_fnn(nilearn_tvt_dir, corr_mat, beh)
        data_ukbb_fnn_sex(nilearn_tvt_dir, corr_mat, beh)
        data_ukbb_brainnetcnn(nilearn_tvt_dir, corr_mat, beh)
        data_ukbb_brainnetcnn_sex(nilearn_tvt_dir, corr_mat, beh)
        data_ukbb_gcnn(nilearn_tvt_dir, corr_mat, beh)
        data_ukbb_gcnn_sex(nilearn_tvt_dir, corr_mat, beh)

        os.makedirs(nilearn_cv_dir, exist_ok=True)
        data_hcp_fnn(nilearn_cv_dir, n_folds, corr_mat, beh)
        data_hcp_brainnetcnn(nilearn_cv_dir, n_folds, corr_mat, beh)
        data_hcp_gcnn(nilearn_cv_dir, n_folds, corr_mat, beh)

        get_gcnn_graph(graph_dir, corr_mat)

        compare_npz_dir(os.path.join(nilearn_cv_dir, 'fnn'), self)
        compare_npz_dir(os.path.join(nilearn_cv_dir, 'brainnetcnn'), self)
        compare_npz_dir(os.path.join(nilearn_cv_dir, 'gcnn'), self)
        compare_npz_dir(nilearn_tvt_dir, self)
        compare_graph_dir(graph_dir, self)

        shutil.rmtree(nilearn_tvt_dir)
        shutil.rmtree(nilearn_cv_dir)
        shutil.rmtree(graph_dir)


class TestUKBBVersionDNN(unittest.TestCase):
    def test_UKBB_fnn(self):
        data_dir = os.path.join(ut_dir, 'results', 'tvt')
        subprocess.call([
            'python3', '../cbig/He2019/CBIG_ukbb_fnn.py', '--path_data',
            data_dir, '--batch_size', '5', '--runs', '1', '--epochs', '5',
            '--gpu', '0'
        ])
        ref = os.path.join(ut_dir, 'results', 'log', 'fnn_pred_1_*.npz')
        result = os.path.join(ut_dir, 'log', 'fnn_pred_1_*.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(result)[0])
        ref = ref['metric']
        res = res['metric']
        self.assertTrue(
            np.isclose(ref, res),
            'UKBB ver. FNN metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'log'))

    def test_UKBB_fnn_sex(self):
        data_dir = os.path.join(ut_dir, 'results', 'tvt')
        subprocess.call([
            'python3', '../cbig/He2019/CBIG_ukbb_fnn_sex.py', '--path_data',
            data_dir, '--batch_size', '5', '--runs', '1', '--epochs', '5',
            '--gpu', '0'
        ])
        ref = os.path.join(ut_dir, 'results', 'log', 'fnn_sex_*.npz')
        result = os.path.join(ut_dir, 'log', 'fnn_sex_*.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(result)[0])
        ref = ref['metric']
        res = res['metric']
        self.assertTrue(
            np.isclose(ref, res),
            'UKBB ver. FNN sex metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'log'))

    def test_UKBB_bcn(self):
        data_dir = os.path.join(ut_dir, 'results', 'tvt')
        subprocess.call([
            'python3', '../cbig/He2019/CBIG_ukbb_brainnetcnn.py',
            '--path_data', data_dir, '--batch_size', '5', '--runs', '1',
            '--epochs', '5', '--dim', '39', '--gpu', '0'
        ])
        ref = os.path.join(ut_dir, 'results', 'log',
                           'brainnetcnn_pred_1_*.npz')
        result = os.path.join(ut_dir, 'log', 'brainnetcnn_pred_1_*.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(result)[0])
        ref = ref['metric']
        res = res['metric']
        self.assertTrue(
            np.isclose(ref, res),
            'UKBB ver. BCN metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'log'))

    def test_UKBB_bcn_sex(self):
        data_dir = os.path.join(ut_dir, 'results', 'tvt')
        subprocess.call([
            'python3', '../cbig/He2019/CBIG_ukbb_brainnetcnn_sex.py',
            '--path_data', data_dir, '--batch_size', '5', '--runs', '1',
            '--epochs', '5', '--dim', '39', '--gpu', '0'
        ])
        ref = os.path.join(ut_dir, 'results', 'log', 'brainnetcnn_sex_*.npz')
        result = os.path.join(ut_dir, 'log', 'brainnetcnn_sex_*.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(result)[0])
        ref = ref['metric']
        res = res['metric']
        self.assertTrue(
            np.isclose(ref, res),
            'UKBB ver. BCN sex metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'log'))

    def test_UKBB_gcnn(self):
        data_dir = os.path.join(ut_dir, 'results', 'tvt')
        graph_dir = os.path.join(ut_dir, 'results', 'graph')
        subprocess.call([
            'python3', '../cbig/He2019/CBIG_ukbb_gcnn.py', '--path_data',
            data_dir, '--runs', '1', '--epochs', '5', '--num_subject', '40',
            '--graph_setup', '40Subject_corr_option_3_param_5',
            '--graph_folder', graph_dir, '--gpu', '0'
        ])
        ref = os.path.join(ut_dir, 'results', 'log', 'gcnn_pred_1_*.npz')
        result = os.path.join(ut_dir, 'log', 'gcnn_pred_1_*.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(result)[0])
        ref = ref['metric']
        res = res['metric']
        print(ref - res)
        self.assertTrue(
            np.isclose(ref, res, rtol=0.1),
            'UKBB ver. GCNN metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'log'))

    def test_UKBB_gcnn_sex(self):
        data_dir = os.path.join(ut_dir, 'results', 'tvt')
        graph_dir = os.path.join(ut_dir, 'results', 'graph')
        subprocess.call([
            'python3', '../cbig/He2019/CBIG_ukbb_gcnn_sex.py', '--path_data',
            data_dir, '--runs', '1', '--epochs', '5', '--num_subject', '40',
            '--graph_setup', '40Subject_corr_option_3_param_5',
            '--graph_folder', graph_dir, '--gpu', '0'
        ])
        ref = os.path.join(ut_dir, 'results', 'log', 'gcnn_sex_*.npz')
        result = os.path.join(ut_dir, 'log', 'gcnn_sex_*.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(result)[0])
        ref = ref['metric']
        res = res['metric']
        print(ref - res)
        self.assertTrue(
            np.isclose(ref, res, rtol=0.1),
            'UKBB ver. GCNN sex metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'log'))


class TestHCPVersionDNN(unittest.TestCase):
    def test_HCP_fnn(self):
        n_folds = config.EXAMPLE_N_FOLDS
        data_dir = os.path.join(ut_dir, 'results', 'cv')
        subprocess.call([
            'python3', '../cbig/He2019/CBIG_hcp_fnn.py', '--path_data',
            data_dir, '--batch_size', '5', '--epochs', '100', '--num_subject',
            '40', '--folds',
            str(n_folds), '--n_measure', '2', '--gpu', '0'
        ])
        ref = os.path.join(ut_dir, 'results', 'log', 'HCP_fnn_*.npz')
        result = os.path.join(ut_dir, 'log', 'HCP_fnn_*.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(result)[0])
        ref = ref['metric']
        res = res['metric']
        self.assertTrue(
            np.isclose(ref, res, rtol=0.5),
            'HCP ver. FNN metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'log'))

    def test_HCP_bcn(self):
        n_folds = config.EXAMPLE_N_FOLDS
        data_dir = os.path.join(ut_dir, 'results', 'cv')
        subprocess.call([
            'python3', '../cbig/He2019/CBIG_hcp_brainnetcnn.py', '--path_data',
            data_dir, '--batch_size', '5', '--epochs', '100', '--num_subject',
            '40', '--folds',
            str(n_folds), '--n_measure', '2', '--gpu', '0'
        ])
        ref = os.path.join(ut_dir, 'results', 'log', 'HCP_brainnetcnn_*.npz')
        result = os.path.join(ut_dir, 'log', 'HCP_brainnetcnn_*.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(result)[0])
        ref = ref['metric']
        res = res['metric']
        self.assertTrue(
            np.isclose(ref, res, rtol=0.5),
            'HCP ver. BCN metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'log'))

    def test_HCP_gcnn(self):
        n_folds = config.EXAMPLE_N_FOLDS
        data_dir = os.path.join(ut_dir, 'results', 'cv')
        graph_dir = os.path.join(ut_dir, 'results', 'graph')
        subprocess.call([
            'python3', '../cbig/He2019/CBIG_hcp_gcnn.py', '--path_data',
            data_dir, '--batch_size', '5', '--epochs', '100', '--num_subject',
            '40', '--folds',
            str(n_folds), '--n_measure', '2', '--graph_setup',
            '40Subject_corr_option_3_param_5', '--graph_folder', graph_dir,
            '--gpu', '0'
        ])
        ref = os.path.join(ut_dir, 'results', 'log', 'HCP_gcnn_*.npz')
        result = os.path.join(ut_dir, 'log', 'HCP_gcnn_*.npz')
        ref = np.load(glob.glob(ref)[0])
        res = np.load(glob.glob(result)[0])
        ref = ref['metric']
        res = res['metric']
        self.assertTrue(
            np.isclose(ref, res, rtol=0.5),
            'HCP ver. GCNN metric ' + str(res) + ' vs. ref ' + str(ref))
        shutil.rmtree(os.path.join(ut_dir, 'log'))


if __name__ == '__main__':
    unittest.main()
