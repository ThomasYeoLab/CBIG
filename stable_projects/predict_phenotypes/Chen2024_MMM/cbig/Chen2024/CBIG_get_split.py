#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import argparse
import numpy as np
import scipy.io as sio
from config import config
from CBIG_misc import read_phes, get_subj_num


def get_split(dataset):
    '''save the split files from KRR classical and meta-matching

    Args:
        dataset (str): the test meta-set dataset name

    Returns:
        None
    '''
    ks = config.KS

    if dataset == 'exp_test': # example data
        input_dir = config.IN_DIR_UT if args.unit_test else config.IN_DIR_EXP
        inter_dir = config.INTER_DIR_UT if args.unit_test else config.INTER_DIR_EXP
        output_dir = config.OUT_DIR_UT if args.unit_test else config.OUT_DIR_EXP
        classical_dir = os.path.join(output_dir, 'output_KRR_classical_' + dataset)
        n_rng = config.N_RNG_EXP
    elif dataset in config.DATASET_NAME['test']:
        classical_dir = os.path.join(config.REP_DIR, 'output_KRR_classical_' + dataset)
        input_dir = config.IN_DIR
        inter_dir = config.INTER_DIR
        n_rng = config.N_RNG
    else:
        raise NameError("Unsupported dataset name")

    os.makedirs(inter_dir, exist_ok=True)
    n_subjects = get_subj_num(input_dir, dataset)

    phes_tes = read_phes(input_dir, dataset)
    split_ind_dict = {}
    for _, phe in enumerate(phes_tes):
        split_ind_dict[phe] = np.zeros((n_rng, len(ks), n_subjects))

    # get split file for k shot
    for _, phe in enumerate(phes_tes):
        for ik, k in enumerate(ks):
            for rng in range(1, n_rng + 1):
                folder = os.path.join(classical_dir, 'output_phe_' + phe, dataset + '_' + phe +
                                      '_k_' + str(k) + '_rng_num_' + str(rng))
                split_mat = os.path.join(folder, dataset + '_subject_split.mat')
                mat = sio.loadmat(split_mat)
                split_ind_dict[phe][rng - 1, ik, :] = np.squeeze(mat['sub_fold']['fold_index'][0][0])
        print(phe)

    for _, phe in enumerate(phes_tes):
        split_ind_dict[phe] = split_ind_dict[phe].astype(bool)

    npz = os.path.join(inter_dir, dataset + '_split_ind_' + str(n_rng) + '.npz')
    np.savez(npz, split_ind_dict=split_ind_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dataset", type=str, default="HCP")

    unit_tests_parser = parser.add_mutually_exclusive_group(required=False)
    unit_tests_parser.add_argument('--unit-test', dest='unit_test', action='store_true')
    unit_tests_parser.add_argument('--not-unit-test', dest='unit_test', action='store_false')
    parser.set_defaults(unit_test=False)

    args = parser.parse_args()
    get_split(args.dataset)
