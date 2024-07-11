#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os


class config:
    BASE_DIR = os.path.join(os.getenv('CBIG_CODE_DIR'),
                            'stable_projects', 'predict_phenotypes', 'Chen2024_MMM')
    REP_DIR = os.path.join(BASE_DIR, 'replication')
    IN_DIR = os.path.join(os.environ['CBIG_REPDATA_DIR'],
                          'stable_projects', 'predict_phenotypes', 'Chen2024_MMM')
    INTER_DIR = os.path.join(REP_DIR, 'output_intermediate')
    OUT_DIR = os.path.join(REP_DIR, 'output')
    MODEL_DIR = os.path.join(REP_DIR, 'models')
    LARGE_DATA_DIR = os.path.join(REP_DIR, 'large_data')

    # meta-data for dataset
    DATASET_NAME = {"extra-large": 'UKBB', 'large': ['ABCD'], 'medium': ['GSP', 'HBN', 'eNKI'],
                    'test': ['HCP', 'HCPA']}  # name of different source datasets
    DATASET_NAME_EXP = {"extra-large": 'exp_train_XL', 'large': ['exp_train_L'], 'medium': ['exp_train_M'],
                    'test': ['exp_test']} # name for example data
    N_RNG = 100 # number of repeat random split
    N_RNG_EXP = 3 # number of repeat random split for example data
    N_CV = 5 # default number of cross-validation folders

    EXP_DIR = os.path.join(BASE_DIR, 'examples')
    IN_DIR_EXP = os.path.join(EXP_DIR, 'exp_input')
    INTER_DIR_EXP = os.path.join(EXP_DIR, 'exp_output', 'output_intermediate')
    OUT_DIR_EXP = os.path.join(EXP_DIR, 'exp_output')
    MODEL_DIR_EXP = os.path.join(EXP_DIR, 'exp_models')

    UT_DIR = os.path.join(BASE_DIR, 'unit_tests')
    IN_DIR_UT = os.path.join(UT_DIR, 'data')
    INTER_DIR_UT = os.path.join(UT_DIR, 'output', 'output_intermediate')
    OUT_DIR_UT = os.path.join(UT_DIR, 'output')
    MODEL_DIR_UT = os.path.join(UT_DIR, 'models')
    REF_OUT_DIR_UT = os.path.join(UT_DIR, 'ref_results')

    BATCH_SIZE = 128
    EPOCHS = 100
    RAMDOM_SEED = 1
    KS = [10, 20, 50, 100, 200]
