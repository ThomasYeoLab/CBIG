#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import glob
import shutil
import zipfile
import argparse
import numpy as np
from joblib import Parallel, delayed

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--in_dir", type=str, required=True)
parser.add_argument("-f", "--fsl_dir", type=str, required=True)
parser.add_argument("-c", "--current_dir", type=str, required=True)
parser.add_argument("-o", "--out_dir", type=str, required=True)
args = parser.parse_args()
print(args)

n_jobs = 23
os.system('source ' + args.fsl_dir + '/etc/fslconf/fsl.sh')
ZIP_TARS = glob.glob(args.in_dir + '/*20227_2_0.zip')
# *20227_2_0.zip is the bulk file directly downloaded from UKBB
np.random.shuffle(ZIP_TARS)  # enforce different order in every instance
print('We have', len(ZIP_TARS), 'zip to run')
os.makedirs(args.out_dir, exist_ok=True)


def conv_fmri2mni(ZIP_TAR):
    '''convert rsfMRI to MNI152 2mm space

    Args:
        ZIP_TAR (str): input 20227 rsfMRI zip file

    Returns:
        None
    '''
    ZIP_OUT_DIR = ZIP_TAR[:-7]
    OUT_MNI_RS = ZIP_TAR.replace(args.in_dir, args.out_dir).replace(
        '20227_2_0.zip', 'MNI_RS.nii.gz')

    if os.path.exists(OUT_MNI_RS):
        return

    try:
        zip_ref = zipfile.ZipFile(ZIP_TAR, 'r')
        zip_ref.extractall(ZIP_OUT_DIR)
        zip_ref.close()

        RS_DIR = os.path.join(ZIP_OUT_DIR, 'fMRI')

        if os.path.exists(RS_DIR + '/unusable'):
            print('unusable!')
        if os.path.exists(RS_DIR + '/rfMRI.ica'):
            print('converting:', ZIP_TAR.split('/')[-1], OUT_MNI_RS)
            msg = args.current_dir + "/applywarp --ref=" + args.current_dir
            msg += "/MNI152_T1_2mm_brain.nii.gz --in=" + RS_DIR
            msg += "/rfMRI.ica/filtered_func_data_clean.nii.gz --out="
            msg += OUT_MNI_RS + " --warp=" + RS_DIR
            msg += "/rfMRI.ica/reg/example_func2standard_warp.nii.gz"
            msg += "--interp=spline"
            os.system(msg)
            # this applywarp is directly copied from fsl/5.0.8/bin
            # I also attached the one I run in this folder
            os.system(
                "%s/fslmaths %s -mas %s/MNI152_T1_2mm_brain_mask_bin.nii.gz %s"
                % (args.current_dir, OUT_MNI_RS, args.current_dir, OUT_MNI_RS))
            # this fslmaths is directly copied from fsl/5.0.8/bin
            # I also attached the one I run of them in this folder
        else:
            print('no rfMRI.ica for', ZIP_TAR[-21:-7])
    except Exception:
        print('Exception!')

    try:
        shutil.rmtree(ZIP_OUT_DIR)
    except Exception:
        pass


if __name__ == '__main__':
    Parallel(n_jobs=n_jobs)(
        delayed(conv_fmri2mni)(ZIP_TAR) for ZIP_TAR in ZIP_TARS[:])
