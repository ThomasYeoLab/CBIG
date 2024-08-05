#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import pandas as pd
import numpy as np
from utils.misc import txt2list, load_pkl
from config import global_config


def read_dataset_pred_results(exp_results_path, models):
    """
    Read dataset prediction results
    """
    exp_results = dict()
    for model in models:
        model_results_path = os.path.join(exp_results_path,
                                          'dataset_pred_' + model + '.txt')
        model_results_vec = np.array(
            txt2list(model_results_path)[:-1]).reshape((-1, 1)).astype(float)
        exp_results[model] = model_results_vec
    return exp_results


def read_glm_results(exp_results_path, models, dataset_pair):
    """
    Read GLM results
    """
    exp_results = dict()
    for model in models:
        model_results_path = os.path.join(
            exp_results_path, 'z_' + model + '_' + dataset_pair + '.csv')
        model_results_df = pd.read_csv(model_results_path)
        exp_results[model] = model_results_df.iloc[:, 1:].values.astype(float)
    return exp_results


def read_manova_results(exp_results_path, models, dataset_pair):
    """
    Read MANOVA results
    """
    exp_results = dict()
    for model in models:
        model_results_path = os.path.join(
            exp_results_path, 'p_' + model + '_' + dataset_pair + '.csv')
        model_results_df = pd.read_csv(model_results_path)
        exp_results[model] = model_results_df.iloc[:, 2].values.astype(float)
    return exp_results


def compare_2array(vec1, vec2, precison=1e-3):
    assert np.allclose(vec1, vec2, atol=precison)


def check_reference_results():
    """
    Check reference results
    """
    results_path = os.path.join(global_config.root_path, 'results')
    ref_pkl_path = os.path.join(global_config.root_path, 'replication',
                                'ref_results', 'CBIG_DeepResBat_ref.pkl')
    ref_pkl = load_pkl(ref_pkl_path)
    dataset_pairs = ['ADNI-AIBL', 'ADNI-MACC']
    experiments = [
        'dataset_pred',
        'glm_DX',
        'manova_DX',
        'glm_MMSE',
        'manova_MMSE',
    ]
    exp_path_dict = {
        'dataset_pred': 'dataset_prediction',
        'glm_DX': 'assoc_glm/DX_demean',
        'manova_DX': 'assoc_manova/DX_demean',
        'glm_MMSE': 'assoc_glm/MMSE_demean',
        'manova_MMSE': 'assoc_manova/MMSE_demean',
    }
    models = [
        'unharm', 'ComBat4cov', 'CovBat4cov', 'cVAE', 'coVAE', 'DeepResBat'
    ]
    for dataset_pair in dataset_pairs:
        # loop over experiments
        for exp in experiments:
            print('..Checking replication for', exp, 'on', dataset_pair)
            exp_results_path = os.path.join(results_path, exp_path_dict[exp],
                                            dataset_pair)
            if exp.startswith('dataset_pred'):
                exp_result = read_dataset_pred_results(exp_results_path,
                                                       models)
            elif exp.startswith('glm'):
                exp_result = read_glm_results(exp_results_path, models,
                                              dataset_pair)
            else:
                exp_result = read_manova_results(exp_results_path, models,
                                                 dataset_pair)
            for model in models:
                print('.....on model:', model)
                compare_2array(ref_pkl[dataset_pair][exp][model],
                               exp_result[model])
    print('You have successfully replicated all results in An2024_DeepResBat!')


if __name__ == '__main__':
    check_reference_results()
