#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import argparse
import subprocess
import numpy as np
from nilearn import datasets
from nilearn import input_data
from nilearn.connectome import ConnectivityMeasure

from cbig.He2019.config import config
from cbig.He2019.CBIG_prepare_data import data_ukbb_fnn, data_ukbb_fnn_sex
from cbig.He2019.CBIG_prepare_data import data_ukbb_brainnetcnn, data_hcp_fnn
from cbig.He2019.CBIG_prepare_data import data_ukbb_gcnn, data_ukbb_gcnn_sex
from cbig.He2019.CBIG_prepare_data import data_hcp_brainnetcnn, data_hcp_gcnn
from cbig.He2019.CBIG_prepare_data import data_ukbb_brainnetcnn_sex
from cbig.He2019.CBIG_prepare_data import get_gcnn_graph


def get_nilearn_adhd_data(n_subjects, nilearn_download_dir):

    # Load the functional datasets
    datasets.get_data_dirs(data_dir=nilearn_download_dir)
    adhd_data = datasets.fetch_adhd(
        n_subjects=n_subjects, data_dir=nilearn_download_dir)
    msdl_data = datasets.fetch_atlas_msdl(data_dir=nilearn_download_dir)
    masker = input_data.NiftiMapsMasker(
        msdl_data.maps,
        resampling_target="data",
        t_r=2.5,
        detrend=True,
        low_pass=.1,
        high_pass=.01,
        memory='nilearn_cache',
        memory_level=1)

    pooled_subjects = []
    adhd_labels = []  # 1 if ADHD, 0 if control
    age = []
    for func_file, confound_file, phenotypic in zip(
            adhd_data.func, adhd_data.confounds, adhd_data.phenotypic):
        time_series = masker.fit_transform(func_file, confounds=confound_file)
        pooled_subjects.append(time_series)
        adhd_labels.append(phenotypic['adhd'])
        age.append(phenotypic['age'])
    correlation_measure = ConnectivityMeasure(kind='correlation')
    corr_mat = correlation_measure.fit_transform(pooled_subjects)
    print('Correlations are stacked in an array of shape {0}'.format(
        corr_mat.shape))
    beh = np.zeros((n_subjects, 2))
    beh[:, 0] = adhd_labels
    beh[:, 1] = age

    return corr_mat, beh


def main(args):
    n_subjects = config.EXAMPLE_N_SUBJECT
    n_folds = config.EXAMPLE_N_FOLDS
    cur_dir = os.path.join(config.CUR_DIR, 'input')
    nilearn_npz = os.path.join(cur_dir, 'nilearn_adhd_fc_beh.npz')

    if args.download or not os.path.isfile(nilearn_npz):
        nilearn_download_dir = os.path.join(config.INTER_DIR, 'nilearn')
        if not os.path.isdir(nilearn_download_dir):
            nilearn_download_dir = cur_dir
        print('Raw nilearn data is downloading at:', nilearn_download_dir)
        print('Processed nilearn data is saving at:', cur_dir)

        corr_mat, beh = get_nilearn_adhd_data(n_subjects, nilearn_download_dir)
        np.savez(nilearn_npz, corr_mat=corr_mat, beh=beh)

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

    gpu = str(args.gpu)
    subprocess.call([
        'python3', '../cbig/He2019/CBIG_ukbb_fnn.py', '--path_data',
        nilearn_tvt_dir, '--batch_size', '5', '--runs', '1', '--epochs', '100',
        '--gpu', gpu
    ])
    subprocess.call([
        'python3', '../cbig/He2019/CBIG_ukbb_fnn_sex.py', '--path_data',
        nilearn_tvt_dir, '--batch_size', '5', '--runs', '1', '--epochs', '100',
        '--gpu', gpu
    ])
    subprocess.call([
        'python3', '../cbig/He2019/CBIG_ukbb_brainnetcnn.py', '--path_data',
        nilearn_tvt_dir, '--batch_size', '5', '--runs', '1', '--epochs', '100',
        '--dim', '39', '--gpu', gpu
    ])
    subprocess.call([
        'python3', '../cbig/He2019/CBIG_ukbb_brainnetcnn_sex.py',
        '--path_data', nilearn_tvt_dir, '--batch_size', '5', '--runs', '1',
        '--epochs', '100', '--dim', '39', '--gpu', gpu
    ])
    subprocess.call([
        'python3', '../cbig/He2019/CBIG_ukbb_gcnn.py', '--path_data',
        nilearn_tvt_dir, '--runs', '1', '--epochs', '1000', '--num_subject',
        '40', '--graph_setup', '40Subject_corr_option_3_param_5',
        '--graph_folder', graph_dir, '--gpu', gpu
    ])
    subprocess.call([
        'python3', '../cbig/He2019/CBIG_ukbb_gcnn_sex.py', '--path_data',
        nilearn_tvt_dir, '--runs', '1', '--epochs', '1000', '--num_subject',
        '40', '--graph_setup', '40Subject_corr_option_3_param_5',
        '--graph_folder', graph_dir, '--gpu', gpu
    ])
    subprocess.call([
        'python3', '../cbig/He2019/CBIG_hcp_fnn.py', '--path_data',
        nilearn_cv_dir, '--batch_size', '5', '--epochs', '100',
        '--num_subject', '40', '--folds',
        str(n_folds), '--n_measure', '2', '--gpu', gpu
    ])
    subprocess.call([
        'python3', '../cbig/He2019/CBIG_hcp_brainnetcnn.py', '--path_data',
        nilearn_cv_dir, '--batch_size', '5', '--epochs', '100',
        '--num_subject', '40', '--folds',
        str(n_folds), '--n_measure', '2', '--gpu', gpu
    ])
    subprocess.call([
        'python3', '../cbig/He2019/CBIG_hcp_gcnn.py', '--path_data',
        nilearn_cv_dir, '--batch_size', '5', '--epochs', '100',
        '--num_subject', '40', '--folds',
        str(n_folds), '--n_measure', '2', '--graph_setup',
        '40Subject_corr_option_3_param_5', '--graph_folder', graph_dir,
        '--gpu', gpu
    ])
    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--gpu', type=int, default=0, help='int, select which gpu to run')
    down_parser = parser.add_mutually_exclusive_group(required=False)
    down_parser.add_argument(
        '--download',
        dest='download',
        action='store_true',
        help='optinal, download the data from nilearn')
    down_parser.add_argument(
        '--no-download',
        dest='download',
        action='store_false',
        help='optinal (default), use saved data from this repository')
    parser.set_defaults(download=False)

    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
