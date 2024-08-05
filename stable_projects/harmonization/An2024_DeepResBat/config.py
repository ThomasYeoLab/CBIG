#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os


class global_config:
    # Path
    root_path = os.path.abspath(os.path.dirname((__file__)))
    data_path = os.path.join(root_path, 'data')
    raw_data_path = os.path.join(data_path, 'raw_data')
    harm_input_path = os.path.join(data_path, 'harm_input')
    columns_path = os.path.join(data_path, 'features', 'columns.txt')
    features_path = os.path.join(data_path, 'features', 'features.txt')
    ROI_features_path = os.path.join(data_path, 'features', 'ROI_features.txt')
    gm_ROI_features = os.path.join(data_path, 'features',
                                   'gm_ROI_features.txt')
    checkpoints_path = os.path.join(root_path, 'checkpoints')
    # ComBat releated
    ComBat_harm_mapper = \
        {'harm_unmatch2match_train_val_ROI': [
            'unmatch2match_train', 'unmatch2match_val'],
            'harm_match2unmatch_train_val_ROI': [
                'match2unmatch_train', 'match2unmatch_val'],
            'harm_unmatch2match_test_ROI': ['unmatch2match_test'],
            'harm_match2unmatch_test_ROI': ['match2unmatch_test'],
            'harm_match2unmatch_train_full_ROI': ['match2unmatch_train_full'],
            'harm_unmatch2match_train_full_ROI': ['unmatch2match_train_full'],
            'harm_match2unmatch_val_full_ROI': ['match2unmatch_val_full'],
            'harm_unmatch2match_val_full_ROI': ['unmatch2match_val_full']}
    # R_HOME_PATH
    env_path = os.environ['CONDA_PREFIX']
    r_home = os.path.join(env_path, 'lib', 'R')
    # GPU
    usingGPU = True
    gpuQ = True
    # HORD releated
    hord_dim = 12
    hord_continous_dim = 5
