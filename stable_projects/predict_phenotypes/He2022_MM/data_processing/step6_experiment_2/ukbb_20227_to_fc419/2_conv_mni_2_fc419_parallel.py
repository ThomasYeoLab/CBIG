#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import glob
import argparse
import numpy as np
from nilearn.input_data import NiftiLabelsMasker
from joblib import Parallel, delayed

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--in_dir", type=str, required=True)
parser.add_argument("-m", "--mask_nii", type=str, required=True)
parser.add_argument("-o", "--out_dir", type=str, required=True)
args = parser.parse_args()
print(args)

masker = NiftiLabelsMasker(labels_img=args.mask_nii, standardize=True)


def corr2_coeff(A, B):
    '''correlation between two matrix

    Args:
        A (ndarray): first matrix for correlation calculation
        B (ndarray): second matrix for correlation calculation

    Returns:
        ndarray, correlation calculated
    '''
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None], ssB[None]))


def compute_419_FC(MNI_TAR):
    '''compute 419 FC from time series in MNI 2mm space

    Args:
        MNI_TAR (str): input time series in MNI 2mm space

    Returns:
        None
    '''

    OUT_FC = os.path.join(args.out_dir, MNI_TAR[-27:-20] + '_alex419_FC.npy')
    ts = masker.fit_transform(MNI_TAR).T
    fc = corr2_coeff(ts, ts)
    np.save(OUT_FC, fc, allow_pickle=False)


if __name__ == '__main__':
    n_jobs = 20
    MNI_TARS = glob.glob(os.path.join(args.in_dir, '*MNI_RS.nii.gz'))
    np.random.shuffle(MNI_TARS)  # enforce different order in every instance
    Parallel(n_jobs=n_jobs)(
        delayed(compute_419_FC)(MNI_TAR) for MNI_TAR in MNI_TARS[:])
