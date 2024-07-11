#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import scipy.io as sio
import argparse
import numpy as np
from config import config


def check_reference_results(args):
    """
    Compare results from replication and refernece
    """

    # compare classical KRR results (for K = 10, 20, 50, 100, 200) on HCP-YA and HCP-Aging
    res_classical_krr = {}
    ref_classical_krr = {}
    ref_classical_krr['HCP'] = {'meta_cor': [0.0299, 0.0475, 0.0788, 0.1118, 0.1562],
                                'meta_cod': [-0.0135, -0.0156, -0.0074, 0.0054, 0.0260]}
    ref_classical_krr['HCPA'] = {'meta_cor': [0.0439, 0.0611, 0.0977, 0.1336, 0.1729],
                                 'meta_cod': [-0.0025, -0.0047, 0.0005, 0.0139, 0.0347]}
    for dataset in ['HCP', 'HCPA']:
        classical_dir = os.path.join(args.rep_dir, 'output_KRR_classical_' + dataset)
        res_classical_krr[dataset] = sio.loadmat(os.path.join(classical_dir, 'final_result', 'krr_classical_res.mat'))
        for metric in ['meta_cor', 'meta_cod']:
            assert np.allclose(
                np.mean(np.mean(res_classical_krr[dataset][metric], axis=0), axis=0),
                ref_classical_krr[dataset][metric],
                atol=1e-3), 'Replication failed for classical KRR results'

    # compare transfer learning results (for K = 10, 20, 50, 100, 200) on HCP-YA and HCP-Aging
    res_transfer_learning = {}
    ref_transfer_learning = {}
    ref_transfer_learning['HCP'] = {'meta_cor': [0.011882, 0.017544, 0.068431, 0.116846, 0.166465],
                                    'meta_cod': [-0.103379, -0.102841, -0.055845, -0.020157, 0.014673]}
    ref_transfer_learning['HCPA'] = {'meta_cor': [0.015951, 0.027117, 0.098746, 0.154228, 0.199809],
                                     'meta_cod': [-0.153972, -0.132331, -0.045206, 0.001887, 0.038298]}
    for dataset in ['HCP', 'HCPA']:
        res_transfer_learning[dataset] = np.load(os.path.join(args.out_dir,
                                                              'transfer_learning_2' + dataset + '_result.npz'))
        for metric in ['meta_cor', 'meta_cod']:
            assert np.allclose(
                np.mean(np.mean(res_transfer_learning[dataset][metric], axis=2), axis=0),
                ref_transfer_learning[dataset][metric],
                atol=1e-3), 'Replication failed for transfer learning'

    # compare meta-matching with stacking results (for K = 10, 20, 50, 100, 200) on HCP-YA and HCP-Aging
    res_MM_stacking = {}
    ref_MM_stacking = {}
    ref_MM_stacking['HCP'] = {'meta_cor': [0.092451, 0.128978, 0.174484, 0.202823, 0.221976],
                              'meta_cod': [0.008966, 0.019855, 0.037773, 0.051762, 0.063005]}
    ref_MM_stacking['HCPA'] = {'meta_cor': [0.123556, 0.156447, 0.191864, 0.213955, 0.232454],
                               'meta_cod': [0.011844, 0.024255, 0.046457, 0.061111, 0.073078]}
    for dataset in ['HCP', 'HCPA']:
        res_MM_stacking[dataset] = np.load(os.path.join(args.out_dir, 'MM_stacking_2' + dataset + '_result.npz'))
        for metric in ['meta_cor', 'meta_cod']:
            assert np.allclose(
                np.mean(np.mean(res_MM_stacking[dataset][metric], axis=2), axis=0),
                ref_MM_stacking[dataset][metric],
                atol=1e-3), 'Replication failed for meta-matching with stacking results'

    # compare meta-matching with dataset stacking results (for K = 10, 20, 50, 100, 200) on HCP-YA and HCP-Aging
    res_dataset_stacking = {}
    ref_dataset_stacking = {}
    ref_dataset_stacking['HCP'] = {'meta_cor': [0.098249, 0.136617, 0.183401, 0.212966, 0.235998],
                                   'meta_cod': [0.010107, 0.020633, 0.039854, 0.055580, 0.069147]}
    ref_dataset_stacking['HCPA'] = {'meta_cor': [0.132396, 0.168016, 0.207661, 0.231073, 0.249811],
                                    'meta_cod': [0.014278, 0.028388, 0.057255, 0.074323, 0.086154]}
    for dataset in ['HCP', 'HCPA']:
        res_dataset_stacking[dataset] = np.load(
            os.path.join(args.out_dir, 'dataset_stacking_2' + dataset + '_result.npz'))
        for metric in ['meta_cor', 'meta_cod']:
            assert np.allclose(
                np.mean(np.mean(res_dataset_stacking[dataset][metric], axis=2), axis=0),
                ref_dataset_stacking[dataset][metric],
                atol=1e-3), 'Replication failed for meta-matching with dataset stacking results'

    # compare multilayer meta-matching results (for K = 10, 20, 50, 100, 200) on HCP-YA and HCP-Aging
    res_multilayer_stacking = {}
    ref_multilayer_stacking = {}
    ref_multilayer_stacking['HCP'] = {'meta_cor': [0.110065, 0.151004, 0.196061, 0.221893, 0.242443],
                                      'meta_cod': [0.012508, 0.024715, 0.044991, 0.059410, 0.071583]}
    ref_multilayer_stacking['HCPA'] = {'meta_cor': [0.134602, 0.170334, 0.209769, 0.233395, 0.253223],
                                       'meta_cod': [0.013996, 0.028640, 0.058358, 0.076118, 0.088461]}
    for dataset in ['HCP', 'HCPA']:
        res_multilayer_stacking[dataset] = np.load(
            os.path.join(args.out_dir, 'multilayer_stacking_2' + dataset + '_result.npz'))
        for metric in ['meta_cor', 'meta_cod']:
            assert np.allclose(
                np.mean(np.mean(res_multilayer_stacking[dataset][metric], axis=2), axis=0),
                ref_multilayer_stacking[dataset][metric],
                atol=1e-3), 'Replication failed for multilayer meta-matching results'

    print('You have successfully replicated all results in Chen2024_MMM!')


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--out_dir', '-o', type=str, default=config.OUT_DIR)
    parser.add_argument('--in_dir', type=str, default=config.IN_DIR)
    parser.add_argument('--inter_dir', type=str, default=config.INTER_DIR)
    parser.add_argument('--model_dir', type=str, default=config.MODEL_DIR)
    parser.add_argument('--rep_dir', type=str, default=config.REP_DIR)

    return parser.parse_args()


if __name__ == '__main__':
    check_reference_results(get_args())
