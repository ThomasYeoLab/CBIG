#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import sys
import random
import shutil
import numpy as np
import pandas as pd

sys.path.append("..")
sys.path.append("../step1_coarse_filter")
sys.path.append("../step3_post_krr_pca")
try:
    from config import config
    from coarse_filter import get_pfc, plot_phe_num_subj
    from post_krr_pca import pca4df, save_df_corrmat_subj_list
except ImportError:
    raise


def pca_train_valid_test():
    '''perform PCA for training and test meta-set

    Args:
        None

    Returns:
        None
    '''

    pca4df(
        pd.read_csv(config.CSV_TRAIN_2ND_FILTER),
        csv_save=config.CSV_TRAIN_3RD_FILTER)
    pca4df(
        pd.read_csv(config.CSV_TEST_2ND_FILTER),
        csv_save=config.CSV_TEST_3RD_FILTER)


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


def split_phes_list(phes, split):
    '''split phenotypes but also keep similar phenotypes in same meta-set

    Args:
        phes (List): list of phenotypes to split
        split (List): target number of phenotypes in each meta-set for split

    Returns:
        List: phenotypes after split
    '''

    phe_unique = np.unique([phe.split('.')[0] for phe in phes])
    phe_unique = [[phe] for phe in phe_unique]
    phe_tmp = [['age-2', '20007-2', '2946-0', '845-2'],
               ['1578-2', '1588-2', '4440-2'], ['20159-0', '23323-2'],
               ['6348-2', '20156-0'], ['398-2', '20133-0'],
               ['31-0', '22022-0']]
    phe_tmp_list = []
    for phe in phe_tmp:
        phe_tmp_list += phe
    for phe in phe_unique:
        if phe[0] not in phe_tmp_list:
            phe_tmp.append(phe)
    phe_unique = phe_tmp
    random.shuffle(phe_unique)
    phe_count = []
    phe_group = []
    for phe in phe_unique:
        cnt = 0
        grp = []
        for tmp in phe:
            cnt += np.sum([i.startswith(tmp) for i in phes])
            grp += [i for i in phes if i.startswith(tmp)]
        phe_count.append(cnt)
        phe_group.append(grp)
    phe_group = np.array(phe_group)

    phe_split = [[], []]
    phe_tmp = phe_group[np.array(phe_count) > 1]
    split_point = int(len(phe_tmp) / 2)
    for phe in phe_tmp[:split_point]:
        for i in phe:
            phe_split[0].append(i)
    for phe in phe_tmp[split_point:]:
        for i in phe:
            phe_split[1].append(i)

    tmp = phe_group[np.array(phe_count) == 1]
    phe_tmp = []
    for phe in tmp:
        for i in phe:
            phe_tmp.append(i)
    need_count = np.array(split) - np.array([len(i) for i in phe_split])
    phe_split[0] += list(phe_tmp[:need_count[0]])
    phe_split[1] += list(phe_tmp[need_count[0]:])
    print(
        'Training meta-set', len(phe_split[0]),
        [config.DICT_PHE_NAME[phe] for phe in phe_split[0]], '\nTest meta-set',
        len(phe_split[1]), [config.DICT_PHE_NAME[phe] for phe in phe_split[1]])

    return phe_split


def save_df_two_meta_set(df1, df2, phes, name, dict_save):
    '''save combined dataframe for meta-matching

    Args:
        df1 (pandas.DataFrame): first data frame
        df2 (pandas.DataFrame): second data frame
        phes (List): list of phenotypes to save
        name (str): name of combined meta-set for this saving
        dict_save (dictionary): store various path for saving

    Returns:
        None
    '''

    df = pd.concat([df1, df2])
    df.sort_values(by=['eid'], inplace=True)
    df.reset_index(drop=True, inplace=True)

    phes.insert(0, 'eid')
    phes.append('25741-2.0')
    if '31-0.0' not in phes:
        phes.append('31-0.0')
    if 'age-2.0' not in phes:
        phes.append('age-2.0')
    df = df[phes]
    save_df_corrmat_subj_list(df, dict_save['csv'], name, dict_save)


def split_phe():
    '''split phenotypes into training and test meta-set

    Args:
        None

    Returns:
        None
    '''

    df_tra = pd.read_csv(config.CSV_TRAIN_3RD_FILTER)
    df_tes = pd.read_csv(config.CSV_TEST_3RD_FILTER)

    del df_tra['Unnamed: 0']
    del df_tes['Unnamed: 0']

    full_phes = list(df_tra.columns.values)
    full_phes.remove('eid')
    full_phes.remove('25741-2.0')
    print('final set of', len(full_phes), 'phenotypes:', full_phes)

    print('\ntrain/test subjects counts:')
    plot_phe_num_subj(
        df_tra, os.path.join(config.DIR_4_OUTPUT, 'num_subj_dist_tra.png'),
        full_phes)
    plot_phe_num_subj(
        df_tes, os.path.join(config.DIR_4_OUTPUT, 'num_subj_dist_tes.png'),
        full_phes)

    seed = config.RAMDOM_SEED
    random.seed(seed)
    np.random.seed(seed)

    phe_tra, phe_tes = split_phes_list(full_phes, [33, 34])

    print('\nsplit phenotypes into:', len(phe_tra), len(phe_tes))
    print(phe_tra, phe_tes)

    save_df_with_phe(df_tra.copy(), list(phe_tra), config.CSV_TRAIN_FINAL,
                     config.PHE_LIST_TRAIN_FINAL)
    save_df_with_phe(df_tes.copy(), list(phe_tes), config.CSV_TEST_FINAL,
                     config.PHE_LIST_TEST_FINAL)

    save_df_two_meta_set(df_tra, df_tes, list(phe_tra), 'train_test',
                         config.DICT_TRAIN_TEST_COMBINE)


def input_process_dnn():
    '''prepare input for the DNN

    Args:
        None

    Returns:
        None
    '''

    seed = config.RAMDOM_SEED
    random.seed(seed)
    np.random.seed(seed)

    npz = np.load(config.DICT_TRAIN_TEST_COMBINE['npz_flat'])
    x_train_raw = npz['pfc_flat'].T
    npz = np.load(config.DICT_SAVE_TEST['npz_flat'])
    tra_csv = pd.read_csv(config.DICT_TRAIN_TEST_COMBINE['csv'])

    fc_tes = npz['pfc_flat'].T
    print('FC train+test shape', x_train_raw.shape)

    # load y for train + validation
    tra_phe = np.genfromtxt(config.PHE_LIST_TRAIN_FINAL, dtype=str)
    y_train_raw = tra_csv[tra_phe].values
    subject_list = tra_csv['eid'].values
    pfc, pfc_flat = get_pfc(subject_list)
    print(np.allclose(pfc_flat.T, x_train_raw))
    print('Y train+test shape', y_train_raw.shape)

    # load x test
    x_test = fc_tes
    print('FC test shape', x_test.shape)

    tes_csv = pd.read_csv(config.CSV_TEST_FINAL)
    tes_phe = np.genfromtxt(config.PHE_LIST_TEST_FINAL, dtype=str)

    y_test_different_set_phe = tes_csv[tes_phe].values
    print('Y test shape', y_test_different_set_phe.shape)
    print(np.sum(~np.isnan(y_test_different_set_phe), 0))

    np.savez(
        config.NPZ_FNN_INPUT,
        x_train_raw=x_train_raw,
        y_train_raw=y_train_raw,
        x_test=x_test,
        y_test_different_set_phe=y_test_different_set_phe,
        tra_phe=tra_phe,
        tes_phe=tes_phe)

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

    shutil.copy(config.NPZ_FNN_INPUT, dir_out)
    shutil.copy(config.CSV_TEST_FINAL, dir_out)
    shutil.copy(config.PHE_LIST_TEST_FINAL, dir_out)
    shutil.copy(config.DICT_SAVE_TEST['mat'], dir_out)
    shutil.copy(config.DICT_SAVE_TEST['subj_list_txt'], dir_out)
    shutil.copy(config.PHE_LIST_TRAIN_FINAL, dir_out)
    shutil.copy(config.DICT_SAVE_TRAIN['subj_list_txt'], dir_out)
    shutil.copy(config.DICT_TRAIN_TEST_COMBINE['csv'], dir_out)
    shutil.copy(config.DICT_TRAIN_TEST_COMBINE['mat'], dir_out)
    shutil.copy(config.DICT_TRAIN_TEST_COMBINE['subj_list_txt'], dir_out)


if __name__ == "__main__":
    os.makedirs(config.DIR_4_OUTPUT, exist_ok=True)
    pca_train_valid_test()
    split_phe()
    input_process_dnn()
    copy_all_files_for_next()
