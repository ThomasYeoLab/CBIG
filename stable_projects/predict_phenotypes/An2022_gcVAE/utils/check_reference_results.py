#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import argparse
import numpy as np
import pandas as pd
from utils.misc import load_pkl, txt2list
from config import global_config


def check_results_args_parser():
    """
    Check result for replication or unittest
    """
    parser = argparse.ArgumentParser(prog='CheckResultsArgs')
    parser.add_argument('--unittest', action='store_true', default=False)

    args, _ = parser.parse_known_args()
    return args


def extract_dataset_pred_results(results_path):
    """
    Extract site prediction results from replication

    Args:
        results_path (str): Path for results
    """
    df_cols = ['Experiment', 'DatasetPair', 'Model', 'DatasetPredAcc']
    dataset_pairs = ['ADNI-AIBL', 'ADNI-MACC']
    rows = []
    # unmatch2match
    result_path = os.path.join(results_path, 'unmatch2match', 'eval_model',
                               'XGBoostPred')
    models = ['unharm', 'ComBat', 'cVAE', 'gcVAE']
    prefix = 'site_pred_'
    for dataset_pair in dataset_pairs:
        for m in models:
            score = float((txt2list(
                os.path.join(result_path, dataset_pair,
                             prefix + m + '.txt'))[-1]).split('_')[0])
            score = round(score, 4)
            row = ['unmatch2match', dataset_pair, m, score]
            rows.append(row)
    # match2match
    result_path = os.path.join(results_path, 'match2unmatch', 'eval_model',
                               'XGBoostPred')
    for dataset_pair in dataset_pairs:
        for m in models:
            score = float((txt2list(
                os.path.join(result_path, dataset_pair,
                             prefix + m + '.txt'))[-1]).split('_')[0])
            score = round(score, 4)
            row = ['match2unmatch', dataset_pair, m, score]
            rows.append(row)
    # sample_size
    result_path = os.path.join(results_path, 'sample_size')
    percs = [
        '10perc', '20perc', '30perc', '40perc', '50perc', '60perc', '70perc',
        '80perc', '90perc'
    ]
    seeds = [
        'seed10', 'seed11', 'seed12', 'seed13', 'seed14', 'seed15', 'seed16',
        'seed17', 'seed18', 'seed19'
    ]
    models = ['ComBat', 'cVAE', 'gcVAE']
    for dataset_pair in dataset_pairs:
        for m in models:
            for perc in percs:
                value = 0
                for seed in seeds:
                    folder_name = perc + '_' + seed
                    score = float((txt2list(
                        os.path.join(result_path, folder_name, 'eval_model',
                                     'XGBoostPred', dataset_pair,
                                     prefix + m + '.txt'))[-1]).split('_')[0])
                    value += score / 10
                value = round(value, 4)
                row = ['sample_size_' + perc, dataset_pair, m, value]
                rows.append(row)
    dataset_pred_df = pd.DataFrame(data=rows, columns=df_cols)
    return dataset_pred_df


def extract_downstream_application_results(results_path, application='DX'):
    """
    Extract downstream application results from replication

    Args:
        results_path (str): Path for results
        application (str, optional): _description_. Defaults to 'DX'.
    """
    if application == 'DX':
        df_cols = ['Experiment', 'DatasetPair', 'Model', 'DiagnosisPredAcc']
    else:
        df_cols = ['Experiment', 'DatasetPair', 'Model', 'MMSEPredMAE']
    rows = []
    dataset_pairs = ['ADNI-AIBL', 'ADNI-MACC']
    # unmatch2match
    result_path = os.path.join(results_path, 'unmatch2match', 'eval_model',
                               'goalDNNPred')
    models = ['ComBat', 'cVAE', 'gcVAE']
    prefix = application + '_pred_result_'
    for dataset_pair in dataset_pairs:
        dataset = dataset_pair.split('-')[1]
        full_models = ['ADNI_unharm', dataset + '_unharm'] + models
        for m in full_models:
            if 'unharm' in m:
                f_name = prefix + m + '.txt'
            else:
                f_name = prefix + dataset + m + '_harm.txt'
            score = float((txt2list(
                os.path.join(result_path, dataset_pair,
                             f_name))[-1]).split('_')[0])
            score = round(score, 4)
            row = ['unmatch2match', dataset_pair, m, score]
            rows.append(row)
    # match2match
    result_path = os.path.join(results_path, 'match2unmatch', 'eval_model',
                               'goalDNNPred')
    for dataset_pair in dataset_pairs:
        dataset = dataset_pair.split('-')[1]
        full_models = ['ADNI_unharm', dataset + '_unharm'] + models
        for m in full_models:
            if 'unharm' in m:
                f_name = prefix + m + '.txt'
            else:
                f_name = prefix + dataset + m + '_harm.txt'
            score = float((txt2list(
                os.path.join(result_path, dataset_pair,
                             f_name))[-1]).split('_')[0])
            score = round(score, 4)
            row = ['match2unmatch', dataset_pair, m, score]
            rows.append(row)
    # sample_size
    result_path = os.path.join(results_path, 'sample_size')
    percs = [
        '10perc', '20perc', '30perc', '40perc', '50perc', '60perc', '70perc',
        '80perc', '90perc'
    ]
    seeds = [
        'seed10', 'seed11', 'seed12', 'seed13', 'seed14', 'seed15', 'seed16',
        'seed17', 'seed18', 'seed19'
    ]
    for dataset_pair in dataset_pairs:
        dataset = dataset_pair.split('-')[1]
        for m in models:
            f_name = prefix + dataset + m + '_harm.txt'
            for perc in percs:
                value = 0
                for seed in seeds:
                    folder_name = perc + '_' + seed
                    score = float((txt2list(
                        os.path.join(result_path, folder_name, 'eval_model',
                                     'goalDNNPred', dataset_pair,
                                     f_name))[-1]).split('_')[0])
                    value += score / 10
                value = round(value, 4)
                row = ['sample_size_' + perc, dataset_pair, m, value]
                rows.append(row)
    pred_df = pd.DataFrame(data=rows, columns=df_cols)
    return pred_df


def check_reference_results(args):
    """
    Compare results from replication and refernece
    Or Compare reuslts from unittest and unittest reference results
    """
    if args.unittest:
        # for unit test
        results_path = os.path.join(global_config.root_path, 'unit_tests',
                                    'results')
        # 1. check dataset prediction
        site_models = ['unharm', 'cVAE', 'gcVAE']
        site_prefix = 'site_pred_'
        site_loger = []
        for site_model in site_models:
            f_path = os.path.join(results_path, 'XGBoost', 'ADNI-AIBL',
                                  site_prefix + site_model + '.txt')
            site_result = float(((txt2list(f_path)[-1]).split('_'))[0])
            site_result = round(site_result, 4)
            site_loger.append(site_result)
        site_loger = np.array(site_loger).astype(float)
        site_pred_ref = np.array(
            [0.9701492537313433, 0.7611940298507462, 0.7761194029850746])
        assert np.allclose(
            site_loger, site_pred_ref,
            atol=1e-3), 'Unit test failed for dataset prediction'
        # 2. check MMSE prediction
        mmse_models = ['AIBL_unharm', 'AIBLcVAE_harm', 'AIBLgcVAE_harm']
        mmse_prefix = 'MMSE_pred_result_'
        mmse_loger = []
        for mmmse_model in mmse_models:
            f_path = os.path.join(results_path, 'goalDNN', 'ADNI-AIBL',
                                  mmse_prefix + mmmse_model + '.txt')
            mmse_result = float(((txt2list(f_path)[-1]).split('_'))[0])
            mmse_result = round(mmse_result, 4)
            mmse_loger.append(mmse_result)
        mmse_loger = np.array(mmse_loger).astype(float)
        mmse_pred_ref = np.array(
            [6.000762594751565, 3.977287407380989, 3.764270954821483])
        assert np.allclose(
            mmse_pred_ref, mmse_loger,
            atol=1e-3), 'Unit test failed for MMSE prediction'
    else:
        # for replications
        ref_pkl = load_pkl(
            os.path.join(global_config.root_path, 'replication', 'ref_results',
                         'CBIG_gcVAE_ref_results.pkl'))
        results_path = os.path.join(global_config.root_path, 'results')
        # 1. compare dataset prediction results
        dataset_pred_df = extract_dataset_pred_results(results_path)
        ref_dataset_pred_df = ref_pkl['dataset_pred']
        assert np.allclose(
            dataset_pred_df.iloc[:, 3].values,
            ref_dataset_pred_df.iloc[:, 3].values,
            atol=1e-3), 'Replication failed for dataset prediction'
        # 2. compare clinical diagnosis prediction results
        dx_pred_df = extract_downstream_application_results(results_path)
        ref_dx_pred_df = ref_pkl['diagnosis_pred']
        assert np.allclose(
            dx_pred_df.iloc[:, 3].values,
            ref_dx_pred_df.iloc[:, 3].values,
            atol=1e-3), 'Replication failed for clinical diagnosis prediction'
        # 3. compare MMSE prediction results
        mmse_pred_df = extract_downstream_application_results(
            results_path, 'MMSE')
        ref_mmse_pred_df = ref_pkl['mmse_pred']
        assert np.allclose(
            mmse_pred_df.iloc[:, 3].values,
            ref_mmse_pred_df.iloc[:, 3].values,
            atol=1e-3), 'Replication failed for MMSE prediction'
        print('You have successfully replicated all results in An2022_gcVAE!')


if __name__ == '__main__':
    check_reference_results(check_results_args_parser())
