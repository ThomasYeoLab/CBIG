#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import pickle
import shutil
import numpy as np
import pandas as pd
from config import global_config


# IO releated functions
def csvs2df(csvs_path, csvs, index_col=None, site=None):
    """
    Read multiple CSV files and return a DataFrame
    """
    for i, csv in enumerate(csvs):
        if i == 0:
            if index_col is None:
                df = pd.read_csv(os.path.join(csvs_path, csv))
            else:
                df = pd.read_csv(os.path.join(csvs_path, csv),
                                 index_col=index_col)
        else:
            if index_col is None:
                temp_df = pd.read_csv(os.path.join(csvs_path, csv))
            else:
                temp_df = pd.read_csv(os.path.join(csvs_path, csv),
                                      index_col=index_col)
            df = df.append(temp_df, ignore_index=True)

    if site is not None:
        return df[df.SITE == site]
    else:
        return df


def df2csv(df, save_path, *header):
    """
    Save a dataframe to a csv file
    """
    if header:
        df.to_csv(save_path, sep=',', index=False, header=header)
    else:
        df.to_csv(save_path, sep=',', index=False)


def save_df(df, save_path, *header):
    """
    Save a dataframe to a csv file
    """
    if header:
        df.to_csv(save_path, sep=',', index=False, header=header)
    else:
        df.to_csv(save_path, sep=',', index=False)


def txt2list(txt_path):
    """
    Read a txt file and return a list
    """
    return [l.strip() for l in open(txt_path)]


def list2txt(list_obj, save_path):
    """
    Save a list to a txt file
    """
    with open(save_path, 'w') as f:
        for item in list_obj:
            f.write("%s\n" % str(item))
    f.close()


def list2csv(list, columns_name, save_path):
    """
    Save a 2D list to csv file
    """
    df = pd.DataFrame(data=list, columns=columns_name)
    df.to_csv(save_path, index=False, sep=',')


def tuple2str(strings_tuple, sperateor='_'):
    """
    Covert all elements in a tuple to a string, joined by seperartor
    """
    return sperateor.join(strings_tuple)


def list2str(strings_list, seperator='_'):
    """
    Covert all elements in a list to a string, joined by seperartor
    """
    return seperator.join(strings_list)


def load_pkl(pkl_path):
    """
    Load a pkl file and return a dictionary
    """
    fobj = open(pkl_path, 'rb')
    data = pickle.load(fobj)
    fobj.close()

    return data


def save_pkl(dict_obj, save_path):
    """
    Save a dictonary to a pkl file
    """
    fobj = open(save_path, 'wb')
    pickle.dump(dict_obj, fobj, protocol=2)
    fobj.close()


def clean_folder(folder_path):
    """
    Delete all files under folder_path
    """
    if os.path.isdir(folder_path):
        shutil.rmtree(folder_path)
    else:
        raise Exception('The passed <folder_path> is not a dircetory')


def clean_nfolds_files(folder_path, files, n_folds=10):
    """
    Delete files of sub folders named as 0, 1, 2, .. under <folder_path>
    """
    for fold in range(n_folds):
        for file in files:
            file_path = os.path.join(folder_path, str(fold), file)
            if os.path.isfile(file_path):
                os.remove(file_path)


def supermakedirs(path, mode=0o777):
    """
    Create folder with given permissions
    """
    if not path or os.path.exists(path):
        return []
    (head, tail) = os.path.split(path)
    res = supermakedirs(head, mode)
    os.mkdir(path)
    os.chmod(path, mode)
    res += [path]
    return res


def create_folder(folder_path, isOverwrite=False):
    """
    Create one folder
    """
    if os.path.isdir(folder_path):
        if isOverwrite:
            clean_folder(folder_path)
            _ = supermakedirs(folder_path, 0o777)
        else:
            pass
    else:
        _ = supermakedirs(folder_path, 0o777)


def create_nfolds_folder(folder_path, nb_folds=10, isOverwrite=True):
    """
    Create sub folders named as 0, 1, 2, .. under <folder_path>
    If exists sub folders, then deleta and recreate by default
    """
    for fold in range(nb_folds):
        sub_folder_path = os.path.join(folder_path, str(fold))
        create_folder(sub_folder_path, isOverwrite)


# Helper functions (Non-IO) used in this project
def cat2levels(values, isDemean=True):
    """
    Convert categorical variables to a multi-level variable
    For eaxmple, if a value vector is [0, 1, 2]
      if isDeman is Falase, should be [[0, 0], [1, 0], [0, 1]];
      if isDemean is True, should bel [[-1/3, -1/3], [2/3, -1/3], [-1/3, 2/3]]
    The first level is whether is MCI or not,
    the second level is whether is AD or not
    """
    values = np.reshape(values, (values.shape[0], ))
    nb_cats = len(np.unique(values))
    multi_level_array = np.zeros((values.shape[0], nb_cats - 1))
    for i in range(1, nb_cats):
        # the first class is implictly expressed
        idx = np.where(values == i)[0]
        multi_level_array[:, i - 1][idx] = 1
    if isDemean:
        mean = np.nanmean(multi_level_array, axis=0)
        return multi_level_array - mean
    else:
        return multi_level_array


def get_bl_df(df):
    """
    Get baseline data of dataframe
    """
    subjects = list(np.unique(df.RID))
    row_index_list = []
    for sub in subjects:
        sub_mask = df.RID == sub
        dates = df[sub_mask]['EXAMDATE'].values
        dates = np.sort(dates)
        bl_date = dates[0]
        bl_mask = (sub_mask) & (df.EXAMDATE == bl_date)
        bl_index = df[bl_mask].index.values[0]
        row_index_list.append(bl_index)
    # only keep baseline data
    bl_df = df.iloc[row_index_list]
    assert bl_df.shape[0] == len(subjects), "bl_df.shape[0] != len(subjects)"
    bl_df.reset_index(inplace=True, drop=True)

    return bl_df


def one_hot(label_vector, nb_classes=2):
    """
    Convert the label_vector to a one-hot vector
    """
    nb_samples = label_vector.shape[0]
    label_vector = label_vector.astype(int)
    one_hot_encoding = np.zeros((nb_samples, nb_classes))
    one_hot_encoding[np.arange(nb_samples), label_vector] = 1

    one_hot_encoding = np.reshape(one_hot_encoding, [-1, nb_classes])

    return one_hot_encoding.astype(float)


def replace_with_harmed_ROI_wrapper(data_path,
                                    harm_out_path,
                                    origin_files,
                                    harmed_files,
                                    nb_folds=10,
                                    sufix=''):
    """
    Wrapper for replacing orginal ROIs with harmonzied ROIs
    """
    assert len(origin_files) == len(harmed_files), \
        'len(origin_files) != len(harmed_files)'
    cols = txt2list(global_config.columns_path)
    for fold in range(nb_folds):
        fold_data_path = os.path.join(data_path, str(fold))
        fold_harm_out_path = os.path.join(harm_out_path, str(fold))
        for i, f in enumerate(origin_files):
            origin_file = f + '.csv'
            harmed_file = 'harm_' + harmed_files[i] + '_ROI' + sufix + '.csv'
            origin_df = pd.read_csv(os.path.join(fold_data_path, origin_file))
            harmed_df = pd.read_csv(
                os.path.join(fold_harm_out_path, harmed_file))
            harmed_origin_df = replace_with_harmed_ROI(origin_df, harmed_df,
                                                       cols)
            harmed_origin_df_save_name = f + sufix + '.csv'
            # save to a csv file
            save_df(
                harmed_origin_df,
                os.path.join(fold_harm_out_path, harmed_origin_df_save_name))


def replace_with_harmed_ROI(origin_df, harmed_df, cols, isTranspose=False):
    """
    Replace orginal ROIs with harmonzied ROIs
    """
    origin_df = origin_df[cols]
    if isTranspose:
        harmed_df = harmed_df.T
    assert origin_df.shape[0] == harmed_df.shape[0], 'Inequal rows(Visits)'
    origin_df.iloc[:, 7:] = harmed_df.iloc[:, :].values

    return origin_df


def add_uniform_noise(vec, eta, seed=0):
    eps = np.random.default_rng(seed=seed).uniform(-eta, eta, size=vec.shape)
    return vec + eps


def add_gaussian_noise(vec, eta, seed=0):
    eps = np.random.default_rng(seed=seed).normal(0, eta, size=vec.shape)
    return vec + eps
