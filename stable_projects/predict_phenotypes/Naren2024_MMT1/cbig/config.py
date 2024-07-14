#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""


class config:
    KS = [10, 20, 50, 100, 200]
    CUTOFF_SIZE = (160, 192, 160)
    CHANNEL_NUMBER = [32, 64, 128, 256, 256, 64]
    OUTPUT_DIM = 67
    OUTPUT_DIM_TRANS = 1
    DROPOUT = 0.1
    IN_C_TRANS = 256
    OUT_C_TRANS = 64
    LAYER_IDX_TRANS = 5
    FLAYER_IDX_TRANS = 6
    WEIGHT_DECAY = 0.02
    MOMENTUM = 0.9
    SCHEDULER_DECREASE = 30
    DATASET = 'ukbb'
    RUNS = 1
    SEED = 1
    C_DIM = 256
    LR = 1e-6
    K_LIMIT = 33
