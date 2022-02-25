#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import sys
import glob
import random
import shutil
import numpy as np
import pandas as pd
import scipy.io as sio

sys.path.append("..")
sys.path.append("../step1_coarse_filter")
sys.path.append("../step3_post_krr_pca")
try:
    from config import config
    from coarse_filter import get_pfc, plot_phe_num_subj
    from post_krr_pca import pca4df, save_df_corrmat_subj_list
except ImportError:
    raise


def load_Schaefer_ROI_ukbb(roi=419):
    '''split subjects into two meta-set and save FC and subject list

    Args:
        roi (int): number of ROIs count

    Returns:
        None
    '''

    csv_diff_roi = config.CSV_DIFF_ROI_UKBB_2ND_FILTER
    if os.path.isfile(csv_diff_roi):
        print('file exists:', csv_diff_roi)
        return

    dir_fc = config.DIR_DIFF_ROI_UKBB[str(roi)]
    subjects = glob.glob(os.path.join(dir_fc, '*.npy'))
    print('We have', len(subjects), 'processed Schaefer ROI FC')

    subjects_real = []
    for i in subjects:
        npz = np.load(i)
        if np.any(np.isnan(npz.flat)):
            print('NAN in', i)
        else:
            subjects_real.append(i)
    subjects = sorted(
        [int(i.split('/')[-1].split('_')[0]) for i in subjects_real])

    print('We have', len(subjects), 'processed Schaefer ROI FC')

    df_phe = pd.read_csv(config.CSV_MAIN_2ND_FILTER)
    print('We have', len(df_phe['eid']), 'subjects in phenotypes CSV')

    subjects = sorted(list(set(subjects) & set(df_phe['eid'].values)))
    print('We have', len(subjects), 'FC matched')

    df = df_phe.loc[df_phe['eid'].isin(subjects)]
    print('df shape after filtered', df.shape)

    save_df_corrmat_subj_list(
        df.copy(),
        csv_diff_roi,
        'schaefer_' + str(roi),
        config.DICT_DIFF_ROI,
        roi=roi)


def pca_diff_roi():
    '''perform PCA for training and test meta-set for experiment 2

    Args:
        None

    Returns:
        None
    '''
    csv_diff_roi = config.CSV_DIFF_ROI_UKBB_3RD_FILTER
    if os.path.isfile(csv_diff_roi):
        print('file exists:', csv_diff_roi)
        return
    pca4df(
        pd.read_csv(config.CSV_DIFF_ROI_UKBB_2ND_FILTER),
        csv_save=csv_diff_roi)


def save_df_with_phe(df, phes, df_name, phes_list_name):
    '''save dataframe with phenotypes

    Args:
        df (pandas.DataFrame): data frame to save
        phes (List): list of phenotypes to save
        df_name (str): path of csv for saving
        phes_list_name (str): path of phenotypes list txt to save

    Returns:
        List: phenotypes after split
    '''

    np.savetxt(phes_list_name, phes, fmt='%s')

    phes.insert(0, 'eid')
    phes.append('25741-2.0')
    if '31-0.0' not in phes:
        phes.append('31-0.0')
    if 'age-2.0' not in phes:
        phes.append('age-2.0')
    df = df[phes]
    df.to_csv(df_name)


def save_final_df_beh_list():
    '''split phenotypes into training and test meta-set

    Args:
        None

    Returns:
        None
    '''

    df_tra = pd.read_csv(config.CSV_DIFF_ROI_UKBB_3RD_FILTER)

    del df_tra['Unnamed: 0']

    full_phes = list(df_tra.columns.values)
    full_phes.remove('eid')
    full_phes.remove('25741-2.0')
    print('final set of', len(full_phes), 'phenotypes:', full_phes)

    print('\ntrain/test subjects counts:')
    plot_phe_num_subj(
        df_tra,
        os.path.join(config.DIR_6_OUTPUT, 'num_subj_dist_across_dataset.png'),
        full_phes)

    save_df_with_phe(df_tra.copy(), list(full_phes), config.CSV_DIFF_ROI_FINAL,
                     config.PHE_LIST_DIFF_ROI_FINAL)


def get_HCP_data():
    '''get all the HCP data needed for experiment 2

    Args:
        None

    Returns:
        Tuple: HCP data subject list, phenotypes list, flatten functional
            connectivity, phenotypes
    '''

    subject_list = np.genfromtxt(config.HCP_SUBJ_LIST, dtype=str)
    print(subject_list.shape, subject_list)

    phe_list = np.genfromtxt(config.HCP_PHE_LIST, dtype=str)
    print(phe_list.shape, phe_list)

    pfc, pfc_flat = get_pfc(subject_list, roi=419, flag_HCP=True)
    print(pfc.shape, pfc_flat.shape)
    sio.savemat(config.HCP_MAT, {'corr_mat': pfc})

    mat = sio.loadmat(
        os.path.join(config.DIR_5_OUTPUT, 'HCP_diff_roi_final_phe.mat'))
    y = mat['y']
    print(y.shape)
    return subject_list, phe_list, pfc_flat, y


def input_process_cross_dataset():
    '''prepare input for the DNN

    Args:
        None

    Returns:
        None
    '''

    seed = config.RAMDOM_SEED
    random.seed(seed)
    np.random.seed(seed)

    npz = np.load(config.DICT_DIFF_ROI['npz_flat'])
    x_train_raw = npz['pfc_flat'].T
    print('FC train shape', x_train_raw.shape)

    # load y for train + validation
    tra_csv = pd.read_csv(config.CSV_DIFF_ROI_FINAL)
    tra_phe = np.genfromtxt(config.PHE_LIST_DIFF_ROI_FINAL, dtype=str)
    y_train_raw = tra_csv[tra_phe].values
    print('Y train shape', y_train_raw.shape)
    print(tra_phe)

    # load x test
    subject_list, phe_list, pfc_flat, y_test = get_HCP_data()
    x_test = pfc_flat.T
    print('FC test shape', x_test.shape)

    y_test_different_set_phe = y_test
    print('Y test shape', y_test_different_set_phe.shape)
    print(np.sum(~np.isnan(y_test_different_set_phe), 0))
    print(phe_list)

    df = pd.DataFrame(y_test, columns=phe_list)
    df['eid'] = subject_list
    df = df[['eid'] + list(phe_list)]
    df.to_csv(config.CSV_HCP_DIFF_ROI_FINAL, index=False)

    np.savez(
        config.NPZ_FNN_INPUT_ACROSS_DATASET,
        x_train_raw=x_train_raw,
        y_train_raw=y_train_raw,
        x_test=x_test,
        y_test_different_set_phe=y_test_different_set_phe,
        tra_phe=tra_phe,
        tes_phe=phe_list)

    return


def copy_all_files_for_next():
    '''copy all the files needed for following experiments

    Args:
        None

    Returns:
        None
    '''

    dir_out = config.DIR_FINAL_OUT
    os.makedirs(dir_out, exist_ok=True)

    shutil.copy(config.HCP_PHE_LIST, dir_out)
    shutil.copy(config.CSV_HCP_DIFF_ROI_FINAL, dir_out)
    shutil.copy(config.HCP_SUBJ_LIST, dir_out)
    shutil.copy(config.HCP_MAT, dir_out)


if __name__ == "__main__":
    os.makedirs(config.DIR_6_OUTPUT, exist_ok=True)
    load_Schaefer_ROI_ukbb()
    pca_diff_roi()
    save_final_df_beh_list()
    input_process_cross_dataset()
    copy_all_files_for_next()
