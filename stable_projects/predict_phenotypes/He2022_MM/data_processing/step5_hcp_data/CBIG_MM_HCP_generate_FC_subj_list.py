#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import sys
import glob

sys.path.append("..")
try:
    from config import config
except ImportError:
    raise


def generate_FC_subj_list():
    '''generate functional connectivity path list file for HCP S1200

    Args:
        None

    Returns:
        None
    '''
    dir_FC = config.DIR_5_FC['419']
    dir_out = config.DIR_5_OUTPUT

    dir_sub = glob.glob(dir_FC + '/*.mat')
    dir_sub = sorted([i.split('_')[-1].split('.')[0] for i in dir_sub])
    print(dir_sub)

    f = open(
        os.path.join(dir_out, 'subject_list_FC_S1200_1094_210120.txt'), "w")

    for i in dir_sub:
        f.write(i)
        f.write('\n')
    f.close()


if __name__ == '__main__':
    generate_FC_subj_list()
