#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import shutil
import pickle
import numpy as np
import pandas as pd
from copy import deepcopy
from config import global_config


class myException(Exception):
    pass


def list2csv(list, columns_name, save_path):
    """
    Save a 2D list to csv file

    Args:
        list (list): List to save
        columns_name (list): Column names of saved csvs
        save_path (str): Path for saving output
    """
    df = pd.DataFrame(data=list, columns=columns_name)
    df.to_csv(save_path, index=False, sep=',')


def list2str(strings_list, seperator='_'):
    """
    Covert all elements in a list to a string, joined by seperartor

    Args:
        strings_list (list): List of strings
        seperator (str, optional): seperator. Defaults to '_'.
    """
    return seperator.join(strings_list)


def txt2list(txt_path):
    """
    Read a txt file and return a list

    Args:
        txt_path (str): Path for test file
    """
    return [element.strip() for element in open(txt_path)]


def list2txt(list_obj, save_path):
    """
    Save a list to a txt file

    Args:
        list_obj (list): List to save
        save_path (str): Path for saving txt
    """
    with open(save_path, 'w') as f:
        for item in list_obj:
            f.write("%s\n" % str(item))
    f.close()


def tuple2str(strings_tuple, sperateor='_'):
    """
    Covert all elements in a tuple to a string, joined by seperartor

    Args:
        strings_tuple (tuple): Tuple to convert
        sperateor (str, optional): sperateor. Defaults to '_'.
    """
    return sperateor.join(strings_tuple)


def csvs2df(csvs_path, csvs, index_col=None, site=None):
    """
    Read multiple CSV files and return a DataFrame

    Args:
        csvs_path (list): List of paths for csv files
        csvs (list): List of csv names
        index_col (_type_, optional): index_col. Defaults to None.
        site (_type_, optional): site. Defaults to None.
    """
    for i, csv in enumerate(csvs):
        if i == 0:
            if index_col is None:
                df = pd.read_csv(os.path.join(csvs_path, csv))
            else:
                df = pd.read_csv(
                    os.path.join(csvs_path, csv), index_col=index_col)
        else:
            if index_col is None:
                temp_df = pd.read_csv(os.path.join(csvs_path, csv))
            else:
                temp_df = pd.read_csv(
                    os.path.join(csvs_path, csv), index_col=index_col)
            df = df.append(temp_df, ignore_index=True)

    if site is not None:
        return df[df.SITE == site]
    else:
        return df


def save_df(df, save_path, *header):
    """
    Save a dataframe to a csv file

    Args:
        df (class DataFrame): DataFrame to save
        save_path (str): Path for saving
        header (str, optional): header
    """
    if header:
        df.to_csv(save_path, sep=',', index=False, header=header)
    else:
        df.to_csv(save_path, sep=',', index=False)


def dfA_substract_dfB(dfA, dfB):
    """
    substract dfB from dfA based on RID & EXAMDATE

    Args:
        dfA (class DataFrame): dfA
        dfB (class DataFrame): dfB
    """
    assert dfA.shape[0] > dfB.shape[0], 'df1.shape[0] <= df2.shape[0]'

    dfA_ = deepcopy(dfA)
    dfB_ = deepcopy(dfB)
    dfA_['key'] = dfA_['RID'] + '-' + dfA_['EXAMDATE']
    dfB_['key'] = dfB_['RID'] + '-' + dfB_['EXAMDATE']

    idx = dfA_['key'].isin(dfB_['key'])
    df = dfA_.drop(dfA_[idx].index)
    df.drop('key', axis=1, inplace=True)

    return df


def load_pkl(pkl_path):
    """
    Load a pkl file and return a dictionary

    Args:
        pkl_path (str): Path for pkl file
    """
    fobj = open(pkl_path, 'rb')
    data = pickle.load(fobj)
    fobj.close()

    return data


def save_pkl(dict_obj, save_path):
    """
    Save a dictonary to a pkl file

    Args:
        dict_obj (dict): dict_obj
        save_path (path): Path for saving dict obj
    """
    fobj = open(save_path, 'wb')
    pickle.dump(dict_obj, fobj, protocol=2)
    fobj.close()


def clean_folder(folder_path):
    """
    Delete all files under folder_path

    Args:
        folder_path (str): Folder path to clean
    """
    if os.path.isdir(folder_path):
        shutil.rmtree(folder_path)
    else:
        raise myException('The passed <folder_path> is not a dircetory')


def clean_nfolds_files(folder_path, files, n_folds=10):
    """
    Delete files of sub folders named as 0, 1, 2, .. under <folder_path>

    Args:
        folder_path (str): Folder path to clean
        files (list): List of files to clean
        n_folds (int, optional): Number of folds. Defaults to 10.
    """
    for fold in range(n_folds):
        for file in files:
            file_path = os.path.join(folder_path, str(fold), file)
            if os.path.isfile(file_path):
                os.remove(file_path)


def create_folder(folder_path, isOverwrite=False):
    """
    Create one folder

    Args:
        folder_path (str): Folder path to clean
        isOverwrite (bool, optional): isOverwrite. Defaults to False.
    """
    if os.path.isdir(folder_path):
        if isOverwrite:
            clean_folder(folder_path)
            os.makedirs(folder_path)
        else:
            pass
    else:
        os.makedirs(folder_path)


def create_nfolds_folder(folder_path, nb_folds=10, isOverwrite=True):
    """
    Create sub folders named as 0, 1, 2, .. under <folder_path>
    If exists sub folders, then deleta and recreate by default

    Args:
        folder_path (str): Folder path to clean
        nb_folds (int, optional): Number of folds. Defaults to 10.
        isOverwrite (bool, optional): isOverwrite. Defaults to False.
    """
    for fold in range(nb_folds):
        sub_folder_path = os.path.join(folder_path, str(fold))
        create_folder(sub_folder_path, isOverwrite)


def one_hot(label_vector, nb_classes=2):
    """
    Convert the label_vector to a one-hot vector

    Args:
        label_vector (ndarray): Vector for labels
        nb_classes (int, optional): Number of classes. Defaults to 2.
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

    Args:
        data_path (str): Path for data
        harm_out_path (str): Path for harmonization output
        origin_files (list): List of names for origin files
        harmed_files (list): List of names for origin files
        nb_folds (int, optional): Number of folds. Defaults to 10.
        sufix (str, optional): Sufix for saving. Defaults to ''.
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

    Args:
        origin_df (class DataFrame): Original dataframe
        harmed_df (class DataFrame): Harmonized dataframe
        cols (list): Columns
        isTranspose (bool, optional): Whether to Transpose. Defaults to False.
    """
    origin_df = origin_df[cols]
    if isTranspose:
        harmed_df = harmed_df.T
    assert origin_df.shape[0] == harmed_df.shape[0], 'Inequal rows(Visits)'
    origin_df.iloc[:, 7:] = harmed_df.iloc[:, :].values

    return origin_df
