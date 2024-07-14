#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import numpy as np
import argparse


def check_results_args_parser():
    '''function to get args from command line and return the args

     Returns:
       argparse.ArgumentParser: args that could be used by other function
    '''

    parser = argparse.ArgumentParser(prog='CheckResultsArgs')
    parser.add_argument('--setup', type=str, default='replication')
    parser.add_argument('--base_dir',
                        type=str,
                        default=os.path.abspath(os.path.join(os.getcwd())))
    parser.add_argument('--res_dir', type=str, default=None)
    parser.add_argument('--ref_dir', type=str, default=None)
    parser.add_argument('--datasetname', type=str, default=None)

    args, _ = parser.parse_known_args()
    return args


def check_reference_results(args):
    """
    Compare results from replication and refernece
    Or Compare reuslts from unittest and unittest reference results
    """
    if args.setup == 'unit_tests':
        # for unit test
        res_dir = os.path.join(args.base_dir, args.setup,
                               'results/unit_tests_test_dataset')
        ref_dir = os.path.join(args.base_dir, args.setup,
                               'ref/unit_tests_test_dataset')
    else:
        # for replications
        res_dir = os.path.join(args.base_dir, args.setup, 'results',
                               args.datasetname)
        ref_dir = os.path.join(args.base_dir, args.setup, 'ref',
                               args.datasetname)

    res_log = []
    ref_log = []

    res = np.load(os.path.join(res_dir,
                               "elasticnet/elasticnet_result_test.npz"),
                  allow_pickle=True)
    res_log.append(
        np.mean(np.mean(np.mean(res.f.meta_cor, axis=2), axis=0), axis=0))
    res_log.append(
        np.mean(np.mean(np.mean(res.f.meta_cod, axis=2), axis=0), axis=0))
    res = np.load(os.path.join(res_dir,
                               "classical_transfer/meta_result_test.npz"),
                  allow_pickle=True)
    res_log.append(
        np.mean(np.mean(np.mean(res.f.meta_cor_tuned, axis=2), axis=0),
                axis=0))
    res_log.append(
        np.mean(np.mean(np.mean(res.f.meta_cod_tuned, axis=2), axis=0),
                axis=0))
    res = np.load(os.path.join(res_dir, "mm_finetune/meta_result_test.npz"),
                  allow_pickle=True)
    res_log.append(
        np.mean(np.mean(np.mean(res.f.meta_cor_tuned, axis=2), axis=0),
                axis=0))
    res_log.append(
        np.mean(np.mean(np.mean(res.f.meta_cod_tuned, axis=2), axis=0),
                axis=0))
    res = np.load(os.path.join(res_dir,
                               "mm_stacking/meta_stacking_result_test.npz"),
                  allow_pickle=True)
    res_log.append(
        np.mean(np.mean(np.mean(res.f.meta_cor, axis=2), axis=0), axis=0))
    res_log.append(
        np.mean(np.mean(np.mean(res.f.meta_cod, axis=2), axis=0), axis=0))

    ref = np.load(os.path.join(ref_dir,
                               "elasticnet/elasticnet_result_test.npz"),
                  allow_pickle=True)
    ref_log.append(
        np.mean(np.mean(np.mean(ref.f.meta_cor, axis=2), axis=0), axis=0))
    ref_log.append(
        np.mean(np.mean(np.mean(ref.f.meta_cod, axis=2), axis=0), axis=0))
    ref = np.load(os.path.join(ref_dir,
                               "classical_transfer/meta_result_test.npz"),
                  allow_pickle=True)
    ref_log.append(
        np.mean(np.mean(np.mean(ref.f.meta_cor_tuned, axis=2), axis=0),
                axis=0))
    ref_log.append(
        np.mean(np.mean(np.mean(ref.f.meta_cod_tuned, axis=2), axis=0),
                axis=0))
    ref = np.load(os.path.join(ref_dir, "mm_finetune/meta_result_test.npz"),
                  allow_pickle=True)
    ref_log.append(
        np.mean(np.mean(np.mean(ref.f.meta_cor_tuned, axis=2), axis=0),
                axis=0))
    ref_log.append(
        np.mean(np.mean(np.mean(ref.f.meta_cod_tuned, axis=2), axis=0),
                axis=0))
    ref = np.load(os.path.join(ref_dir,
                               "mm_stacking/meta_stacking_result_test.npz"),
                  allow_pickle=True)
    ref_log.append(
        np.mean(np.mean(np.mean(ref.f.meta_cor, axis=2), axis=0), axis=0))
    ref_log.append(
        np.mean(np.mean(np.mean(ref.f.meta_cod, axis=2), axis=0), axis=0))

    res_log, ref_log = np.array(res_log), np.array(ref_log)
    assert np.allclose(res_log, ref_log, atol=1e-3), 'Failed for prediction'

    print('You have successfully replicated all results in Naren2024_MMT1!')


if __name__ == "__main__":
    check_reference_results(check_results_args_parser())
