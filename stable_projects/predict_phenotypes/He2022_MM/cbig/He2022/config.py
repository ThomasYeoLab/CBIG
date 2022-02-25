#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os


class config:
    BASE_DIR = os.path.join(
        os.getenv('CBIG_CODE_DIR'),
        'stable_projects/predict_phenotypes/He2022_MM')
    REP_DIR = os.path.join(BASE_DIR, 'replication')
    IN_DIR = os.path.join(REP_DIR, 'input_dnn')
    INTER_DIR = os.path.join(REP_DIR, 'output_intermediate')
    OUT_DIR = os.path.join(REP_DIR, 'output_dnn')

    LARGE_DATA_DIR = os.path.join('/home/the/storage/data/metalearning')

    BATCH_SIZE = 128
    EPOCHS = 1000
    RUNS = 1
    RAMDOM_SEED = 1
    KS = [10, 20, 50, 100, 200]
