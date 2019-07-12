#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import numpy as np


class config:
    BASE_DIR = '../../../../../../data/fmri_predict_behavior'
    CUR_DIR = os.getcwd()
    INTER_DIR = os.path.join(BASE_DIR, 'He2019_data')
    GRAPH_FOLDER = os.path.join(INTER_DIR, 'graph')
    RAMDOM_SEED = 42
    OUT_PATH = 'log'

    # Config for HCP
    HCP_CORR_MAT = 'FC_subject_953.mat'
    HCP_SUBJECT_LIST = 'He2019_hcp_953_split.mat'
    HCP_ORIG_DIR = os.path.join(BASE_DIR, 'original_data_953')
    HCP_INTER_DIR = os.path.join(INTER_DIR, 'HCP')
    HCP_MEASURE_SETS = ['Cognitive', 'Personality_Task', 'Social_Emotion']
    HCP_NUM_FOLD = 20
    HCP_NUM_SUBJECT = 953
    HCP_N_DIMENSION = 419
    HCP_BATCH_SIZE = 128
    HCP_MEASURE_SETS_NUM = [13, 22, 23]
    HCP_N_MEASURE = int(np.sum(HCP_MEASURE_SETS_NUM))

    # Config for UKBB
    UKBB_CORR_MAT = 'ukbb_ht_180205_FC_55.mat'
    UKBB_SUBJECT_LIST = 'ukbb_subject_split.mat'
    UKBB_ORIG_DIR = os.path.join(BASE_DIR, 'original_data_ukbb_8868')
    UKBB_INTER_DIR = os.path.join(INTER_DIR, 'UKBB')
    UKBB_MEASURE_SETS = ['1802_8868']
    UKBB_NUM_SUBJECT = 8868
    UKBB_RUNS = 5
    UKBB_BATCH_SIZE = 128
    UKBB_EPOCHS = 200
    UKBB_EPOCHS_GCNN = 2000
    UKBB_N_DIM = 55

    # Config for example
    EXAMPLE_N_SUBJECT = 40
    EXAMPLE_N_FOLDS = 4
