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


def generate_TC_subj_list():
    '''generate time series path list file for HCP S1200

    Args:
        None

    Returns:
        None
    '''
    dir_base = config.DIR_HCP_BASE
    dir_sub = glob.glob(dir_base + '/*')
    dir_sub = [os.path.join(i, 'MNINonLinear', 'Results') for i in dir_sub]
    scans_name = [
        'rfMRI_REST1_LR', 'rfMRI_REST1_RL', 'rfMRI_REST2_LR', 'rfMRI_REST2_RL'
    ]

    f = open(config.SUBJ_LIST_HCP_TC, "w")
    rnt = []
    for i in dir_sub:
        tmp = []
        for scan in scans_name:
            scan_file = os.path.join(
                i, scan, scan + '_Atlas_MSMAll_hp2000_clean.dtseries.nii')
            if os.path.isfile(scan_file):
                tmp.append(scan_file)
        if len(tmp) > 0:
            for j in tmp:
                f.write(j + ' ')
            f.write('\n')
            rnt.append(tmp)
    f.close()


if __name__ == '__main__':
    os.makedirs(config.DIR_5_OUTPUT, exist_ok=True)
    generate_TC_subj_list()
