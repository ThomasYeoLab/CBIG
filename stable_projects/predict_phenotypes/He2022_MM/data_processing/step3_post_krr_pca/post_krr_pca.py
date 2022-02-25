#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import sys
import random
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib import rcParams
from statsmodels.multivariate.pca import PCA
from sklearn.model_selection import train_test_split

sys.path.append("..")
sys.path.append("../step1_coarse_filter")
try:
    from config import config
    from coarse_filter import get_pfc, plot_phe_num_subj
except ImportError:
    raise

plt.switch_backend('Agg')
rcParams['font.family'] = 'Arial'
pd.set_option('display.max_columns', None)


def print_phenotype_info(phe_name, phe_info, phe_data):
    '''print the info of selected phenotypes, info includes data-field id,
    name, different level of category

    Args:
        phe_name (List): list of selected phenotypes
        phe_info (pandas.DataFrame): UK Biobank showcase information DataFrame
        phe_data (pandas.DataFrame): data of each phenotypes

    Returns:
        None
    '''

    print(len(phe_name), 'numbers of phenotypes have been selected')
    phe_unique = []
    for i in phe_name:
        temp = i.split('-')[0]
        if temp == 'age':
            print('age is also included.')
        else:
            phe_unique.append(int(i.split('-')[0]))
    indexes = np.unique(phe_unique, return_index=True)[1]
    phe_unique = [phe_unique[index] for index in sorted(indexes)]
    print(len(phe_unique), 'numbers of UNIQUE phenotypes have been selected')

    full_list = []
    phe_info_new = phe_info[phe_info['data_id'].isin(phe_unique)]
    category = phe_info_new['level1'].unique()
    for cat in category:
        phe_info_temp = phe_info_new[phe_info_new['level1'] == cat]
        print(cat)
        level2 = phe_info_temp['level2'].unique()
        for l2 in level2:
            phe_info_l2 = phe_info_temp[phe_info_temp['level2'] == l2]
            print('   ', l2)
            level3 = phe_info_l2['level3'].unique()
            for l3 in level3:
                phe_info_l3 = phe_info_l2[phe_info_l2['level3'] == l3]
                if l3 == '0.0':
                    for ind, row in phe_info_l3.iterrows():
                        tmp = [
                            phe for phe in phe_name
                            if phe.startswith(str(row['data_id']) + '-')
                        ]
                        tmp_subj_cnt = [
                            np.sum(~np.isnan(phe_data[phe].values))
                            for phe in tmp
                        ]
                        full_list = full_list + tmp
                        print('       ', row['data_id'], row['data_field'],
                              tmp, tmp_subj_cnt)
                else:
                    print('       ', l3)
                    for ind, row in phe_info_l3.iterrows():
                        tmp = [
                            phe for phe in phe_name
                            if phe.startswith(str(row['data_id']) + '-')
                        ]
                        tmp_subj_cnt = [
                            np.sum(~np.isnan(phe_data[phe].values))
                            for phe in tmp
                        ]
                        full_list = full_list + tmp
                        print('           ', row['data_id'], row['data_field'],
                              tmp, tmp_subj_cnt)
    print(len(full_list), full_list)


def save_phe_list(phes, phes_removed):
    print(len(phes), 'numbers of phenotypes have been selected')
    phe_unique = []
    for i in phes:
        temp = i.split('-')[0]
        if temp == 'age':
            print('age is also included.')
        else:
            phe_unique.append(int(i.split('-')[0]))
    indexes = np.unique(phe_unique, return_index=True)[1]
    phe_unique = [phe_unique[index] for index in sorted(indexes)]
    print(len(phe_unique), 'numbers of UNIQUE phenotypes have been selected')
    print(sorted(phe_unique), 'age')

    print(len(phes_removed), 'numbers of phenotypes have been removed')
    phe_removed_unique = []
    for i in phes_removed:
        temp = i.split('-')[0]
        if temp == 'age':
            print('age is also included.')
        else:
            phe_removed = int(i.split('-')[0])
            if phe_removed not in phe_unique:
                phe_removed_unique.append(int(i.split('-')[0]))
    indexes = np.unique(phe_removed_unique, return_index=True)[1]
    phe_removed_unique = [
        phe_removed_unique[index] for index in sorted(indexes)
    ]
    print(
        len(phe_removed_unique),
        'numbers of UNIQUE phenotypes have been removed')
    print(sorted(phe_removed_unique))


def check_phe_select_krr_output():
    '''check output from KRR result for phenotypes

    Args:
        None

    Returns:
        None
    '''

    # get data
    npz = np.load(config.NPZ_PHE_SELECT_RES, allow_pickle=True)
    phe_corr = npz['final_corr_avg']
    phe_name = npz['final_name']
    phe_corr_df = pd.DataFrame({'id': phe_name, 'corr': phe_corr})
    phe_corr_df['id'] = [i.split('-')[0] for i in phe_name]
    phe_info = pd.read_csv(config.CSV_UKBB_CAT_INFO)
    phe_data = pd.read_csv(config.CSV_COARSE_FILTER)
    phe_info = phe_info[[
        'data_id', 'data_field', 'Strata', 'Array', 'level1', 'level2',
        'level3', 'level4', 'level5', 'level6'
    ]]
    phe_info.columns = [
        'data_id', 'data_field', 'level1', 'level2', 'level3', 'level4',
        'level5', 'level6', 'level7', 'level8'
    ]

    # use 0.1 as threshold and display all phenotypes in their own category
    # can refer to section 4.4 of meta-matching paper
    threshold = 0.1
    phe_ind = phe_corr > threshold
    phes = phe_name[phe_ind]
    phes_removed = phe_name[phe_corr <= threshold]
    phe_corr = phe_corr[phe_ind]

    # print information
    save_phe_list(phes, phes_removed)
    print_phenotype_info(phes, phe_info, phe_data)

    # save csv with new set of phenotype after KRR selection
    if 'eid' not in phes:
        print('eid not in phes')
        phes = np.append(phes, 'eid')
    if 'age-2.0' not in phes:
        print('age-2.0 not in phes')
        phes = np.append(phes, 'age-2.0')
    if '25741-2.0' not in phes:
        print('25741-2.0 not in phes')
        phes = np.append(phes, '25741-2.0')
    if '31-0.0' not in phes:
        print('31-0.0 not in phes')
        phes = np.append(phes, '31-0.0')
    print(phes)
    dat_csv = pd.read_csv(config.CSV_1000_FOR_PHE_SELECT)
    filtered_csv = dat_csv.loc[:, phes]
    filtered_csv = filtered_csv.sort_values(by=['eid'])
    filtered_csv.to_csv(config.CSV_1000_2ND_FILTER)

    dat_csv = pd.read_csv(config.CSV_REMAIN)
    filtered_csv = dat_csv.loc[:, phes]
    filtered_csv = filtered_csv.sort_values(by=['eid'])
    filtered_csv.to_csv(config.CSV_MAIN_2ND_FILTER)

    # plot hist of number of subject for each phenotype
    df_phe = pd.read_csv(config.CSV_1000_FOR_PHE_SELECT)
    plot_phe_num_subj(df_phe,
                      os.path.join(config.DIR_3_OUTPUT, 'num_subj_dist.png'),
                      phes)


def save_df_corrmat_subj_list(df, csv, name, dict_save, roi=55):
    '''save csv, FC, subject list for meta-set dataframe

    Args:
        df (pandas.DataFrame): the meta-set data frame to be saved
        csv (str): path of csv for saving
        name (str): name of meta-set for this saving
        dict_save (dictionary): store various path for saving

    Returns:
        None
    '''

    df.sort_index(inplace=True)
    if not np.array_equal(df['eid'], np.sort(df['eid'])):
        raise Exception(name + 'df eid is not sorted')
    if 'Unnamed: 0' in df.columns.values:
        del df['Unnamed: 0']

    plot_phe_num_subj(
        df, os.path.join(config.DIR_3_OUTPUT,
                         'num_subj_dist_' + name + '.png'))
    df.reset_index(drop=True, inplace=True)
    df.to_csv(csv)

    subject_list = df['eid'].values
    pfc, pfc_flat = get_pfc(subject_list, roi=roi)

    np.savez(dict_save['npz'], pfc=pfc)
    np.savez(dict_save['npz_flat'], pfc_flat=pfc_flat)
    if roi == 55:
        sio.savemat(dict_save['mat'], dict([('corr_mat', pfc)]))
        sio.savemat(dict_save['mat_flat'], dict([('pfc_flat', pfc_flat)]))
    np.savez(dict_save['subj_list'], subject_list=subject_list)
    np.savetxt(dict_save['subj_list_txt'], subject_list, fmt='%s')


def split_subject_meta():
    '''split subjects into two meta-set and save FC and subject list

    Args:
        None

    Returns:
        None
    '''

    num_val_tes = 10000
    # refer to section 4.6 of meta-matching paper

    csv_train = config.CSV_TRAIN_2ND_FILTER
    csv_test = config.CSV_TEST_2ND_FILTER
    if os.path.isfile(csv_train) and os.path.isfile(csv_test):
        print('file exists:', csv_train, 'and', csv_test)
        return

    df_phe = pd.read_csv(config.CSV_MAIN_2ND_FILTER)
    seed = config.RAMDOM_SEED
    random.seed(seed)
    np.random.seed(seed)
    df_train, df_test = train_test_split(
        df_phe, test_size=num_val_tes, random_state=seed)

    save_df_corrmat_subj_list(df_train.copy(), csv_train, 'train',
                              config.DICT_SAVE_TRAIN)
    save_df_corrmat_subj_list(df_test.copy(), csv_test, 'test',
                              config.DICT_SAVE_TEST)


def combine_column(df, phes, ncomp=6, need_PCA=False):
    '''Perform colume combine for dataframe on certain phenotypes (columns)

    Args:
        df (pandas.DataFrame): the data frame to perform colume combine on
        phes (list): list of phenotypes to perform colume combine
        ncomp (int, optional): number of component for PCA, default is 6
        need_PCA (bool, optional): perform PCA if true, otherwise just average

    Returns:
        Tuple: dataframe with combined columes, list of combined columes name
    '''

    # get value to perform colume combine on from input dataframe
    value = np.zeros((len(df.index), len(phes)))
    for i, phe in enumerate(phes):
        value[:, i] = df[phe].values
        df.drop(columns=phe, inplace=True)
    ind_real = np.sum(~np.isnan(value), axis=1).astype(bool)
    ind_all_real = np.sum(~np.isnan(value), axis=1) == len(phes)
    ind_ncomp = np.sum(~np.isnan(value), axis=1) >= ncomp
    print(ind_ncomp.sum(), ind_real.sum(), ind_all_real.sum())

    value = value[ind_ncomp, :]

    # combine columes with PCA or average
    if need_PCA:
        x = PCA(value, ncomp=ncomp, standardize=True, missing='fill-em')
        variance = x.rsquare
        i = 0
        new_phes = []
        while i < ncomp:
            new_phe = phes[0].split('.')[0] + '.' + str(999 + i)
            new_phes.append(new_phe)
            processed_value = np.empty((len(df.index)))
            processed_value[:] = np.nan
            processed_value[ind_ncomp] = x.factors[:, i]
            df.loc[:, new_phe] = processed_value
            i += 1
        print('pca to', len(new_phes), 'phes:',
              phes[0].split('.')[0] + '.999 to', new_phe,
              'with variance explained', variance[i])

    else:
        processed_value = np.nanmean(value, axis=1)
        new_phe = phes[0].split('.')[0] + '.999'
        print(phes, 'combined to', new_phe)
        df[new_phe] = processed_value

    return df, new_phes


def pca4df(df, csv_save=None, csv_pca_phe=None):
    '''Perform PCA for dataframe

    Args:
        df (pandas.DataFrame): the data frame for PCA
        csv_save (str, optional): path of csv for saving
        csv_pca_phe (str, optional): path of PCA phenotype csv for next KRR
            filter

    Returns:
        None
    '''

    print('number of subjects:', len(df.index))

    phes_pca = []

    phes = [
        '20156-0.0', '20157-0.0', '20149-0.1', '20149-0.7', '20149-0.13',
        '20149-0.14', '20149-0.17', '20149-0.25', '20155-0.4', '20155-0.7',
        '20155-0.15', '20155-0.16', '20155-0.21'
    ]  # Trail making
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=6)
    phes_pca = phes_pca + phes_new

    phes = [
        '20159-0.0', '20195-0.0', '20200-0.22', '20200-0.33', '20229-0.33',
        '20230-0.3', '20230-0.6', '20230-0.13', '20230-0.15', '20230-0.16',
        '20230-0.18', '20230-0.26'
    ]  # Symbol digit substitution
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=6)
    phes_pca = phes_pca + phes_new

    phes = ['22022-0.0', '22023-0.0']  # Genotyping
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=2)
    phes_pca = phes_pca + phes_new

    phes = ['22003-0.0', '22009-0.1', '22009-0.11']  # Genotyping
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=2)
    phes_pca = phes_pca + phes_new

    phes = ['20007-2.0', '40008-0.0']  # cancer
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=1)
    phes_pca = phes_pca + phes_new

    phes = ['30512-1.0', '30502-1.0', '30522-1.0', '30532-1.0']  # Urine assays
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=1)
    phes_pca = phes_pca + phes_new

    phes = [
        '30620-0.0', '30620-1.0', '30630-0.0', '30630-1.0', '30650-0.0',
        '30650-1.0', '30700-0.0', '30700-1.0', '30720-0.0', '30720-1.0',
        '30730-0.0', '30740-1.0', '30750-0.0', '30760-0.0', '30760-1.0',
        '30770-0.0', '30770-1.0', '30790-1.0', '30800-0.0', '30830-0.0',
        '30830-1.0', '30850-0.0', '30850-1.0', '30840-1.0', '30870-0.0',
        '30880-0.0', '30880-1.0', '30670-0.0', '30030-0.0', '30030-1.0',
        '30030-2.0', '30020-0.0', '30020-1.0', '30020-2.0', '30300-0.0',
        '30300-1.0', '30300-2.0', '30290-1.0', '30280-1.0', '30180-2.0',
        '30050-2.0', '30040-1.0', '30040-2.0', '30270-2.0', '30080-0.0',
        '30080-1.0', '30080-2.0', '30090-0.0', '30090-2.0', '30010-0.0',
        '30010-1.0', '30010-2.0', '30250-0.0', '30250-1.0', '30250-2.0',
        '30240-1.0', '30162-1.0', '30222-1.0', '30152-1.0', '30212-1.0',
        '30032-1.0', '30022-1.0', '30302-2.0', '30292-2.0', '30282-2.0',
        '30122-1.0', '30182-1.0', '30052-1.0', '30062-1.0', '30042-1.0',
        '30102-1.0', '30262-1.0', '30262-2.0', '30272-2.0', '30132-1.0',
        '30192-1.0', '30142-1.0', '30202-1.0', '30082-1.0', '30012-1.0',
        '30072-1.0', '30252-1.0', '30252-2.0', '30242-1.0', '30242-2.0',
        '30002-1.0'
    ]  # Blood assays
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=6)
    phes_pca = phes_pca + phes_new

    phes = ['26414-0.0', '26410-0.0']  # Multiple Deprivation
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=2)
    phes_pca = phes_pca + phes_new

    phes = [
        '21671-2.0', '21663-2.0', '21651-2.0', '21621-2.0', '21664-2.0',
        '3-2.0', '21631-2.0'
    ]  # Process durations
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=5)
    phes_pca = phes_pca + phes_new

    phes = [
        '6348-2.0', '6772-2.17', '6772-2.22', '6773-2.2', '6773-2.5',
        '6773-2.6', '6773-2.14', '6773-2.24'
    ]  # Trail making
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=6)
    phes_pca = phes_pca + phes_new

    phes = ['21004-2.0', '6382-2.0']  # Tower rearranging
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=1)
    phes_pca = phes_pca + phes_new

    phes = ['23323-2.0', '23324-2.0']  # Symbol digit substitution
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=1)
    phes_pca = phes_pca + phes_new

    phes = ['20023-2.0', '404-2.0', '404-2.1', '404-2.5',
            '404-2.11']  # Reaction time
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=3)
    phes_pca = phes_pca + phes_new

    phes = ['4286-2.0', '4288-2.0', '4289-2.0']  # Prospective memory
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=2)
    phes_pca = phes_pca + phes_new

    phes = ['4250-2.7', '4253-2.2', '4253-2.7', '4255-2.6',
            '4256-2.6']  # Numeric memory
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=3)
    phes_pca = phes_pca + phes_new

    phes = ['6373-2.0', '6374-2.0', '6333-2.0',
            '6333-2.7']  # Matrix pattern completion
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=4)
    phes_pca = phes_pca + phes_new

    phes = [
        '20009-2.1', '20009-2.2', '20009-2.3', '20008-2.4', '135-2.0',
        '137-2.0'
    ]  # non-cancer illness
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=5)
    phes_pca = phes_pca + phes_new

    phes = [
        '20075-2.0', '130-1.0', '129-0.0', '24508-0.0', '22702-0.3',
        '22704-0.0', '22704-0.2'
    ]  # geometry
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=3)
    phes_pca = phes_pca + phes_new

    phes = [
        '3062-2.0', '3062-2.1', '3062-2.2', '20151-0.0', '3063-2.0',
        '3063-2.1', '3063-2.2', '20150-0.0', '20153-0.0', '3064-2.0',
        '3064-2.1', '3064-2.2'
    ]  # Spirometry
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=4)
    phes_pca = phes_pca + phes_new

    phes = ['46-2.0', '47-2.0']  # Hand grip strength
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=2)
    phes_pca = phes_pca + phes_new

    phes = [
        '12336-2.0', '12340-2.0', '5986-1.30', '5986-1.31', '5986-1.32',
        '5986-1.33', '5986-1.34', '5986-1.35', '5986-1.36', '5986-1.37',
        '5986-1.38', '5986-1.41', '5986-1.42', '5986-1.43', '5986-1.44',
        '5986-1.48', '5986-1.49', '5986-1.50', '5986-1.51', '5986-1.52',
        '5986-1.56', '5986-1.57', '5986-1.58', '5986-1.73', '5984-1.4',
        '5984-1.5', '5984-1.6', '5984-1.7', '5984-1.8', '5984-1.9',
        '5984-1.10', '5984-1.11', '5984-1.12', '5984-1.13', '5984-1.14',
        '5984-1.15', '5984-1.16', '5984-1.17', '5984-1.18', '5984-1.19',
        '5984-1.20', '5984-1.21', '5984-1.22', '5984-1.23', '5984-1.24',
        '5984-1.25', '5984-1.26', '5984-1.27', '5984-1.28', '5984-1.29',
        '5984-1.30', '5984-1.31', '5984-1.32', '5984-1.33', '5984-1.34',
        '5984-1.35', '5984-1.36', '5984-1.37', '5984-1.38', '5984-1.39',
        '5984-1.40', '5984-1.41', '5984-1.42', '5984-1.43', '5984-1.45',
        '5984-1.46', '5984-1.47', '5984-1.48', '5984-1.49', '5984-1.50',
        '5984-1.51', '5984-1.52', '5984-1.53', '5984-1.54', '5984-1.55',
        '5984-1.56', '5984-1.57', '5984-1.58', '5984-1.59', '5984-1.60',
        '5984-1.61', '5984-1.62', '5984-1.63', '5984-1.66', '5984-1.67',
        '5984-1.68', '5984-1.69', '5984-1.73', '5984-1.74', '5983-1.8',
        '5983-1.25', '5983-1.26', '5983-1.27', '5983-1.28', '5983-1.29',
        '5983-1.30', '5983-1.31', '5983-1.33', '5983-1.34', '5983-1.35',
        '5983-1.36', '5983-1.37', '5983-1.38', '5983-1.39', '5983-1.40',
        '5983-1.41', '5983-1.42', '5983-1.43', '5983-1.44', '5983-1.45',
        '5983-1.46', '5983-1.47', '5983-1.48', '5983-1.49', '5983-1.50',
        '5983-1.51', '5983-1.52', '5983-1.53', '5983-1.54', '5983-1.68',
        '5983-1.69', '5983-1.70', '5983-1.71', '5983-1.72', '5983-1.76',
        '6033-1.0', '6032-1.0'
    ]  # ECG
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=6)
    phes_pca = phes_pca + phes_new

    phes = [
        '22672-2.0', '22675-2.0', '22678-2.0', '22681-2.0', '22671-2.0',
        '22674-2.0', '22677-2.0', '22680-2.0', '22670-2.0', '22673-2.0',
        '22676-2.0', '22679-2.0'
    ]  # Carotid ultrasound
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=6)
    phes_pca = phes_pca + phes_new

    phes = [
        '3143-0.0', '3144-0.0', '3147-0.0', '3148-0.0', '78-0.0', '3085-0.0',
        '77-0.0', '4100-2.0', '4119-2.0', '4101-2.0', '4120-2.0', '4104-2.0',
        '4123-2.0', '4105-2.0', '4124-2.0', '4106-2.0', '4125-2.0', '4146-2.0',
        '4144-2.0', '4145-2.0', '4138-2.0', '4143-2.0'
    ]  # Bone-densitometry of heel
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=3)
    phes_pca = phes_pca + phes_new

    phes = [
        '12144-2.0', '12143-2.0', '48-2.0', '21002-2.0', '21001-2.0', '49-2.0',
        '50-2.0', '20015-2.0', '23098-2.0', '23104-2.0', '23113-2.0',
        '23118-2.0', '23114-2.0', '23123-2.0', '23119-2.0', '23124-2.0',
        '23120-2.0', '23121-2.0', '23125-2.0', '23126-2.0', '23122-2.0',
        '23127-2.0', '23128-2.0', '23129-2.0', '23130-2.0', '23105-2.0',
        '23099-2.0', '23100-2.0', '23101-2.0', '23102-2.0', '23115-2.0',
        '23111-2.0', '23116-2.0', '23112-2.0', '23117-2.0', '23106-2.0',
        '23110-2.0', '23109-2.0', '23108-2.0', '23107-2.0'
    ]  # Anthropometry
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=4)
    phes_pca = phes_pca + phes_new

    phes = [
        '4194-2.0', '5115-1.0', '5100-1.0', '5306-1.0', '5109-1.0', '5109-1.1',
        '5109-1.2', '5157-1.0', '5157-1.1', '5157-1.2', '5114-1.0', '5162-1.0',
        '5162-1.1', '5106-1.1', '5101-1.0', '5089-1.0', '5257-1.0', '5262-1.0',
        '5263-1.0', '4079-2.0', '4079-2.1', '94-2.0', '94-2.1', '95-2.0',
        '95-2.1', '102-2.0', '102-2.1', '4080-2.0', '4080-2.1', '93-2.0',
        '93-2.1'
    ]  # blood and eye measures
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=6)
    phes_pca = phes_pca + phes_new

    phes = ['2946-0.0', '2946-2.0', '1845-2.0']  # father/mother age
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=2)
    phes_pca = phes_pca + phes_new

    phes = ['20162-2.0', '20161-2.0']  # smoking
    df, phes_new = combine_column(df, phes, need_PCA=True, ncomp=1)
    phes_pca = phes_pca + phes_new

    df_copy = df.copy(deep=True)

    phe_ind = [
        '4440-2.0', '1578-2.0', '1588-2.0', '864-2.0', '1090-2.0', '1070-2.0',
        '1160-2.0', '845-2.0', '767-2.0', '777-2.0', '709-2.0', '20127-0.0',
        '4230-2.3', '20016-2.0', '398-2.3', '31-0.0', '20133-0.2', 'age-2.0'
    ]  # individual phenotypes that not going to be processed with PCA

    if csv_save:
        npz = np.load(config.NPZ_PHE_PCA_SELECT_RES, allow_pickle=True)
        tmp_corr_avg = npz['final_corr_avg']
        tmp_name = npz['final_name']
        tmp = tmp_name[tmp_corr_avg > 0.1]

        final_phes = phe_ind + list(tmp)
        print(len(final_phes), final_phes)
        save_phes = final_phes.copy()
        save_phes.insert(0, 'eid')
        save_phes.append('25741-2.0')
        df = df[save_phes]

        df.to_csv(csv_save)

    if csv_pca_phe:
        print('we have ', len(phes_pca), 'pca generated phenotypes')
        np.savetxt(config.PHE_LIST_PCA, phes_pca, fmt='%s')

        phes_pca.insert(0, 'eid')
        phes_pca.append('25741-2.0')
        phes_pca.append('31-0.0')
        phes_pca.append('age-2.0')
        df_copy = df_copy[phes_pca]
        df_copy.to_csv(csv_pca_phe)


def pca_1000_subject():
    '''perform PCA for 1000 subject

    Args:
        None

    Returns:
        None
    '''

    pca4df(
        pd.read_csv(config.CSV_1000_2ND_FILTER),
        csv_pca_phe=config.CSV_1000_PCA_PHE)


if __name__ == "__main__":
    os.makedirs(config.DIR_3_OUTPUT, exist_ok=True)
    check_phe_select_krr_output()
    split_subject_meta()
    pca_1000_subject()
