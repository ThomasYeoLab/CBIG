#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import sys
import random
import datetime
import numpy as np
import pandas as pd
import scipy.io as sio
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
plt.switch_backend('Agg')

sys.path.append("..")
try:
    from config import config
except ImportError:
    raise


def filter_subject_by_FC():
    '''filter the raw CSV file by only include subjects with FC ('25751-2.0')

    Args:
        None

    Returns:
        None
    '''

    raw_csv_name = config.CSV_RAW

    csv_base_name = config.CSV_BASE
    if not os.path.isfile(csv_base_name):
        if 'raw_csv' not in locals():
            raw_csv = pd.read_csv(raw_csv_name)

        tmp = raw_csv['25751-2.0'].dropna()
        print(tmp.shape[0], 'numbers of subjects with FC by UKBB')
        csv_base = raw_csv.loc[tmp.index.values, :]
        csv_base.to_csv(csv_base_name)
    else:
        print('file exists:', csv_base_name)


def convert_date2float(data):
    '''convert dataframe's date to float

    Args:
        data (pandas.DataFrame): df that to convert the date to float

    Returns:
        pandas.DataFrame: df with date converted to float
    '''
    data[~data.isnull()] = pd.to_datetime(
        data.dropna()).astype(int).astype(float)
    return data.astype(float)


def unique_entry(x):
    '''check number of unique entry

    Args:
        x (list): entry list

    Returns:
        int: number of unique entry
    '''
    return len(np.unique([tmp.split('-')[0] for tmp in x]))


def coarse_filter():
    '''perform coarse filter to phenotypes based on various criterion, like:
    remove category phenotypes except sex
    remove brain MRI phenotypes
    remove phenotypes have less than 2000 participants
    remove phenotypes with same value for more than 80% participants
    (details in section 4.4 in the meta-matching paper)

    Args:
        None

    Returns:
        None
    '''

    if os.path.isfile(config.CSV_COARSE_FILTER):
        print('file exists:', config.CSV_COARSE_FILTER)
        return

    # load data
    dat_csv = pd.read_csv(config.CSV_BASE)
    ori_csv = dat_csv.copy()
    inf_csv = pd.read_csv(config.CSV_UKBB_CAT_INFO)

    # check entry in the UKBioBank showcase
    # '20235-0.0' is deleted here due to not in showcase and all nan value
    entry_v0 = dat_csv.columns.values[4:]
    print('number of entry before selection:', len(entry_v0), 'unique:',
          unique_entry(entry_v0))
    inf_list = inf_csv['data_id'].values
    entry_v1 = []
    for entry in entry_v0:
        phenotype = int(entry.split('-')[0])
        if phenotype not in inf_list:
            print(entry, 'is not in the UKBioBank showcase webstie')
        else:
            entry_v1.append(entry)
    print('number of entry check ukbb showcase:', len(entry_v1), 'unique:',
          unique_entry(entry_v1))

    # Removed non-continuous and non-integer data fields (date and time
    # converted to float)
    # Removed Brain MRI phenotypes (category ID 100)
    entry_v2 = []
    for entry in entry_v1:
        phenotype = int(entry.split('-')[0])
        row = inf_csv.loc[inf_csv['data_id'] == phenotype]
        phe_type = row['Value Type'].values[0]
        phe_brain_mri = row['level1'].values[0]
        if phe_brain_mri == 'Brain MRI':
            # print(entry, row['data_field'].values[0])
            continue
        if phe_type.startswith('Date') or phe_type.startswith('Time'):
            dat_csv[entry] = convert_date2float(dat_csv[entry].copy())
            entry_v2.append(entry)
        if phe_type.startswith('Continuous') or phe_type.startswith('Integer'):
            entry_v2.append(entry)
    print('number of entry after keep Continuous or Integer value:',
          len(entry_v2), 'unique:', unique_entry(entry_v2))

    # only keep first imaging visit, except early visit have more subjects
    entry_v2_5 = entry_v2
    entry_v2_5 = np.array(entry_v2_5)
    phe_list = np.unique([int(i.split('-')[0]) for i in entry_v2_5])
    for i in phe_list:

        # Removed first repeat imaging visit (instance 3)
        temp_list = [j for j in entry_v2_5 if j.startswith(str(i) + '-')]
        for j in temp_list:
            if j.startswith(str(i) + '-3.'):
                entry_v2_5 = np.delete(entry_v2_5,
                                       np.where(entry_v2_5 == j)[0][0])
        temp_list = [j for j in temp_list if not j.startswith(str(i) + '-3.')]

        # Removed first two instances (instance 0 and 1) if first imaging visit
        # (instance 2) exists and first imaging visit (instance 2) participants
        # were more than double of participants from instance 0 or 1
        temp_list_2 = np.sort(np.unique([j.split('.')[0] for j in temp_list]))
        if len(temp_list_2) in [2, 3] and str(i) + '-2' == temp_list_2[-1]:
            n_subj = [0, 0, 0]
            for j in temp_list:
                tmp = dat_csv[j].dropna().values.shape[0]
                tmp_i = int(j.split('.')[0].split('-')[1])
                if n_subj[tmp_i] < tmp:
                    n_subj[tmp_i] = tmp
            print(i, n_subj)
            if n_subj[2] * 2 >= n_subj[1] and n_subj[2] * 2 >= n_subj[0]:
                for j in temp_list:
                    if not j.startswith(str(i) + '-2'):
                        entry_v2_5 = np.delete(entry_v2_5,
                                               np.where(entry_v2_5 == j)[0][0])

        # Removed first instance (instance 0) if only the first two instances
        # (instance 0 and 1) exist, and instance 1 participants were more than
        # double of participants from instance 0
        if len(temp_list_2) == 2 and str(i) + '-1' == temp_list_2[-1]:
            n_subj = [0, 0]
            for j in temp_list:
                tmp = dat_csv[j].dropna().values.shape[0]
                tmp_i = int(j.split('.')[0].split('-')[1])
                if n_subj[tmp_i] < tmp:
                    n_subj[tmp_i] = tmp
            print(i, n_subj)
            if n_subj[1] * 2 >= n_subj[0]:
                for j in temp_list:
                    if not j.startswith(str(i) + '-1'):
                        entry_v2_5 = np.delete(entry_v2_5,
                                               np.where(entry_v2_5 == j)[0][0])
    print('number of entry after remove previous visit:', len(entry_v2_5),
          'unique:', unique_entry(entry_v2_5))

    # Removed bulk item
    # Removed phenotypes for which less than 2000 participants had RSFC data
    # Removed behaviors with the same value for more than 80% of participants
    entry_v3 = []
    for entry in entry_v2_5:
        phe_temp = entry.replace('-', '_').replace('.', '_')
        phe_value = dat_csv[entry].dropna().values
        phe_type = dat_csv[entry].dtype
        _, phe_count = np.unique(phe_value, return_counts=True)
        if phe_temp in phe_value:
            print(entry)
        elif np.logical_not(np.all(np.isreal(phe_value))):
            print(entry, 'not real')
        elif phe_type == object:
            print(entry)
        elif phe_value.shape[0] < 2000:
            print(entry, phe_value.shape[0])
        elif np.max(phe_count) / np.sum(phe_count) > 0.8:
            print(entry, phe_count)
        else:
            entry_v3.append(entry)
    print('after delete bulk item:', len(entry_v3), 'unique:',
          unique_entry(entry_v3))

    entry_v3.append('31-0.0')  # add back sex

    # calculate detailed age
    entry_v3.append('age-2.0')
    entry_v3.remove('34-0.0')
    entry_v3.remove('53-2.0')
    entry_v3.remove('21022-0.0')

    np.savetxt(
        os.path.join(config.DIR_1_OUTPUT, 'phe_continuous_new.txt'),
        entry_v3,
        fmt='%s')

    entry_v3.remove('age-2.0')
    entry_v3.append('eid')
    entry_v3.append('25741-2.0')
    final_csv = dat_csv.loc[:, entry_v3]
    temp = final_csv.loc[:, 'eid'].astype(int)
    final_csv.loc[:, 'eid'] = temp

    age = np.zeros((final_csv.shape[0]))
    for i in range(final_csv.shape[0]):
        year = int(ori_csv['34-0.0'].values[i])
        month = int(ori_csv['52-0.0'].values[i])
        scan_date = ori_csv['53-2.0'].values[i]
        scan_date = datetime.datetime.strptime(scan_date, '%Y-%m-%d').date()
        birth_date = datetime.datetime.strptime(
            str(year) + str(month), '%Y%m').date()
        age[i] = (scan_date - birth_date).days / 365

    final_csv['age-2.0'] = pd.Series(age, index=final_csv.index)

    # Save out filtered CSV
    for i in final_csv.columns.values[1:]:
        print(i, np.sum(~np.isnan(final_csv[i].values)))

    final_csv = final_csv.sort_values(by=['eid'])

    final_phes = list(final_csv.columns.values)
    final_phes.remove('eid')
    final_phes.remove('25741-2.0')
    print('phenotypes count after first filter', len(final_phes), 'unique:',
          unique_entry(final_phes))

    final_csv.to_csv(config.CSV_COARSE_FILTER)


def plot_phe_num_subj(df_phe, save_name, phes=None):
    '''plot hist of dataframe subjects count for various phenotypes

    Args:
        df_phe (pandas.DataFrame): phenotype DataFrame
        save_name (str): save place of png files
        phes (List): phenotypes to plot

    Returns:
        None
    '''

    subj_cnt = []
    if phes is None:
        phes = df_phe.columns.values
    for phe in phes:
        cnt = np.sum(~np.isnan(df_phe[phe].values))
        subj_cnt.append(cnt)
        # print(phe, np.sum(~np.isnan(df_phe[phe].values)))
    fig, ax = plt.subplots()
    ax.hist(subj_cnt, 50)
    ax.set(xlabel='num of subjects', ylabel='num of phenotype')
    plt.tight_layout()
    fig.savefig(save_name)
    plt.close('all')

    phe_count = []
    for phe in phes:
        phe_count.append(df_phe[phe].dropna().values.shape[0])
    phe_count = np.sort(phe_count)
    print(save_name, df_phe.shape, phe_count)


def split_for_phe_select():
    '''split 1000 subjects for next step phenotype selection

    Args:
        None

    Returns:
        None
    '''

    csv_1000 = config.CSV_1000_FOR_PHE_SELECT
    csv_remain = config.CSV_REMAIN
    if os.path.isfile(csv_1000) and os.path.isfile(csv_remain):
        print('file exists:', csv_1000, 'and', csv_remain)
        return

    df_phe = pd.read_csv(config.CSV_COARSE_FILTER)
    seed = config.RAMDOM_SEED
    random.seed(seed)
    np.random.seed(seed)
    df_train, df_test = train_test_split(
        df_phe, test_size=1000, random_state=seed)

    df_train.sort_index(inplace=True)
    df_test.sort_index(inplace=True)

    if not np.array_equal(df_test['eid'], np.sort(df_test['eid'])):
        raise Exception('df test eid is not sorted')
    if not np.array_equal(df_train['eid'], np.sort(df_train['eid'])):
        raise Exception('df train eid is not sorted')

    del df_train['Unnamed: 0']
    del df_test['Unnamed: 0']
    plot_phe_num_subj(
        df_train, os.path.join(config.DIR_1_OUTPUT,
                               'num_subj_dist_remain.png'))

    phe_count = []
    for phe in df_test.columns.values:
        phe_count.append(df_test[phe].dropna().values.shape[0])
    phe_count = np.sort(phe_count)
    print(phe_count)

    df_test.reset_index(drop=True, inplace=True)
    df_test.to_csv(csv_1000)
    df_train.reset_index(drop=True, inplace=True)
    df_train.to_csv(csv_remain)

    plot_phe_num_subj(
        pd.read_csv(csv_1000),
        os.path.join(config.DIR_1_OUTPUT, 'num_subj_dist_1000.png'))


def get_pfc(subject_list, roi=55, flag_HCP=False):
    '''convert txt FC file downloaded form UK Biobank to ndarray

    Args:
        subject_list (ndarray): list of subjects eid
        roi (int): number of roi for FC
        flag_HCP (optional, bool): whether get data for HCP dataset

    Returns:
        Tuple: ndarray for partial FC, one is in NxN matrix, one is flatten
    '''

    index = np.tril(np.ones(roi), k=-1) == 1
    if roi == 55:
        dir_pfc = config.DIR_PFC
    else:
        if flag_HCP:
            dir_pfc = config.DIR_5_FC[str(roi)]
        else:
            dir_pfc = config.DIR_DIFF_ROI_UKBB[str(roi)]
    pfc_list = []
    pfc_list_flat = []
    for i in subject_list:
        if roi == 55:
            temp = np.zeros((roi, roi))
            temp_flat = np.loadtxt(
                os.path.join(dir_pfc,
                             str(i) + '_25753_2_0.txt'))
            temp[index] = temp_flat
            temp = temp + temp.T
        else:
            if flag_HCP:
                temp = sio.loadmat(
                    os.path.join(
                        dir_pfc,
                        'FC_' + str(roi) + '_ROIs_' + str(i) + '.mat'))
                temp = temp['corr_mat']
            else:
                temp = np.load(
                    os.path.join(dir_pfc,
                                 str(i) + '_alex' + str(roi) + '_FC.npy'))
            # print(temp.shape)
            temp_flat = temp[index]
            # print(temp_flat.shape)
        pfc_list.append(temp)
        pfc_list_flat.append(temp_flat)
    pfc = np.stack(pfc_list, axis=2)
    pfc_flat = np.stack(pfc_list_flat, axis=1)

    return pfc, pfc_flat


def create_pfc_corr_mat():
    '''create mat file for partial FC for 1000 subjects
    These 1000 subjects are going to be used only for phenotypes selection

    Args:
        None

    Returns:
        None
    '''
    # get partial FC and subject list

    if os.path.isfile(config.MAT_PFC_1000_FLAT) and os.path.isfile(
            config.PHE_LIST_COARSE_FILTER):
        print('file exists:', config.MAT_PFC_1000_FLAT)
        print('file exists:', config.PHE_LIST_COARSE_FILTER)
        return

    csv = pd.read_csv(config.CSV_1000_FOR_PHE_SELECT)
    subject_list = csv['eid'].values
    pfc, pfc_flat = get_pfc(subject_list)

    np.savez(config.NPZ_PFC_1000, pfc=pfc)
    sio.savemat(config.MAT_PFC_1000, dict([('corr_mat', pfc)]))
    np.savez(config.NPZ_PFC_1000_FLAT, pfc_flat=pfc_flat)
    sio.savemat(config.MAT_PFC_1000_FLAT, dict([('pfc_flat', pfc_flat)]))
    np.savez(config.SUBJ_LIST_1000, subject_list=subject_list)
    np.savetxt(config.SUBJ_LIST_1000_TXT, subject_list, fmt='%s')

    phe_list = csv.columns.values
    phe_list = np.delete(phe_list, -2)
    phe_list = np.delete(phe_list, -2)
    phe_list = np.delete(phe_list, 0)
    print('we have ', len(phe_list), 'phenotypes')
    np.savetxt(config.PHE_LIST_COARSE_FILTER, phe_list, fmt='%s')


if __name__ == "__main__":
    os.makedirs(config.DIR_1_OUTPUT, exist_ok=True)
    filter_subject_by_FC()
    coarse_filter()
    split_for_phe_select()
    create_pfc_corr_mat()
