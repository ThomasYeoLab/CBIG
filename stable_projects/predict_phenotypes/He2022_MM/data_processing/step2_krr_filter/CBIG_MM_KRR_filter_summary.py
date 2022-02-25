#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import sys
import argparse
import numpy as np
import scipy.io as sio
sys.path.append("..")
try:
    from config import config
except ImportError:
    raise


def check_running(rng_num, base_data_dir, phe_list):
    """check whether all phenotypes have KRR results

    Args:
        rng_num (int): number of rng
        base_data_dir (str): directory of KRR result to look at.
        phe_list (str): path of subject list

    Returns:
        int: number of phenotypes have KRR results
    """
    result_folders = []
    for phe in phe_list:
        result_folders.append(
            os.path.join(base_data_dir, 'output_phe_' + phe,
                         'ukbb_' + phe + '_rng_num_' + str(rng_num) + '/'))
    result_folders.sort()

    count = 0
    print('RNG num:', str(rng_num),
          ', phenotypes measures have not finished correctly:')
    for result_folder in result_folders:
        if os.path.isfile(result_folder + 'final_result.mat'):
            count += 1
        else:
            print(result_folder.split('_')[5])

    return count


def load_corr(rng_num, base_data_dir, phe_list):
    """load results from output of the KRR phe filter of specific rng_num

    Args:
        rng_num (int): number of rng
        base_data_dir (str): directory of KRR result to look at.
        phe_list (str): path of subject list

    Returns:
        Tuple: array of phenotypes KRR prediction correlation for rng_num
            and array of phenotypes name
    """

    result_files = []
    for phe in phe_list:
        folder = os.path.join(base_data_dir, 'output_phe_' + phe,
                              'ukbb_' + phe + '_rng_num_' + str(rng_num))
        result_files.append(os.path.join(folder, 'final_result.mat'))

    temp_corr = np.zeros((len(phe_list)))
    temp_name = []

    for i, result_file in enumerate(result_files):
        temp = sio.loadmat(result_file)
        optimal_corr = temp['optimal_acc']
        optimal_corr = optimal_corr[0][0]
        temp_corr[i] = optimal_corr
        temp_name.append(result_file.split('_')[-5])

    return temp_corr, np.array(temp_name)


def summarize_output(phe_list, base_data_dir, output_npz):
    """summarize the output of the KRR phe filter

    Args:
        phe_list (str): path of subject list
        base_data_dir (str): directory of KRR result to look at.
        output_npz (str): path of output npz

    Returns:
        None
    """

    flag_check = 0
    flag_finished = 1
    num_of_rng_num = 100
    phe_list = np.genfromtxt(phe_list, dtype='str')

    # remove 20131-0.2 because it only have 1 value
    temp = np.argwhere(phe_list == '20131-0.2')
    phe_list = np.delete(phe_list, temp)

    num_of_phe = len(phe_list)

    if flag_check:  # flag to check whether KRR have result
        for rng in range(1, num_of_rng_num + 1):
            temp_num = check_running(rng, base_data_dir, phe_list)
            if temp_num != num_of_phe:
                print('rng_num' + str(rng) + ' only has ' + str(temp_num) +
                      ' phenotypes')
                flag_finished = 0

    if flag_finished:  # flag state whether all KRR have finished running
        final_corr = np.zeros((num_of_phe, num_of_rng_num))
        final_name = np.zeros((num_of_phe, num_of_rng_num)).astype(object)

        for rng in range(1, num_of_rng_num + 1):
            final_corr[:, rng - 1], final_name[:, rng - 1] = load_corr(
                rng, base_data_dir, phe_list)
        for i in range(1, num_of_phe):
            temp = final_name[i, 0]
            for j in range(0, num_of_rng_num):
                if temp != final_name[i, j]:
                    print(temp + ' not match with' + final_name[i, j])

        final_corr_avg = np.nanmean(final_corr[:, 0:], axis=1)

        final_name = final_name[:, 0]
        np.savez(
            output_npz,
            final_corr_avg=final_corr_avg,
            final_corr=final_corr,
            final_name=final_name)

        return


def get_args():
    '''function to get args from command line and return the args
    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()
    pca_parser = parser.add_mutually_exclusive_group(required=False)
    pca_parser.add_argument('--pca', dest='pca', action='store_true')
    pca_parser.add_argument('--not-pca', dest='pca', action='store_false')
    parser.set_defaults(pca=False)

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    if args.pca:
        summarize_output(config.PHE_LIST_PCA, config.DIR_2_OUTPUT_PCA,
                         config.NPZ_PHE_PCA_SELECT_RES)
    else:
        summarize_output(config.PHE_LIST_COARSE_FILTER, config.DIR_2_OUTPUT,
                         config.NPZ_PHE_SELECT_RES)
