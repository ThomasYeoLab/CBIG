#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import numpy as np
import pandas as pd
from munkres import Munkres
from copy import deepcopy
from utils.misc import load_pkl, list2csv, create_folder


def hangarian_matching_wrapper(args, bin, step4_output_path):
    """
    Wrapper function for hangarian matching

    Args:
        args (tuple): Parameters
        bin (int): Bin
        step4_output_path (str): Path for saving step4 output
    """
    # get step3_dist_path
    step3_dist_path = os.path.join(
        args.checkpoint_path, args.matching_pair,
        'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
        'round_' + str(args.round), 'step3_distance',
        'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
        + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty))
    threshold_path = os.path.join(step4_output_path,
                                  'Threshold_' + str(args.threshold))
    create_folder(threshold_path)
    # step3-1: generate squared distiance matrix
    for tp in range(args.threshold, 0, -1):
        # matching tps from high to low
        dist_pkl_name = str(args.matching_pair) + '_' + \
                        'AGE' + str(args.age_penalty) + '_' + \
                        'SEX' + str(args.sex_penalty) + '_' + \
                        'DX' + str(args.dx_penalty) + '_' + \
                        'MMSE' + str(args.mmse_penalty) + '_' + \
                        'TP' + str(tp) + '.pkl'
        dist = load_pkl(os.path.join(step3_dist_path, dist_pkl_name))
        # generate squared cost matrix and index looking up table
        gen_squared_cost_matrix(dist, threshold_path, tp)
    # step3-2: Run hangrain matching algorithm
    hangarian(args, step3_dist_path, threshold_path, args.threshold)
    # step3-3: merge matched tables
    merge_matched_tables(threshold_path, args.threshold)


def gen_squared_cost_matrix(dist, output_path, tp, name1='MACC', name2='ADNI'):
    """
    Generate squared cost matrix for each tp and save under output_path

    Args:
        dist (dict): Dictionary for distance between dataset1 and dataset2
        output_path (str): Path for saving output
        tp (int): Number of time points to begin matching
        name1 (str, optional): Name for dataset1. Defaults to 'MACC'.
        name2 (str, optional): Name for dataset2. Defaults to 'ADNI'.
    """
    nb_subs_dataset1 = len(dist.keys())
    nb_subs_dataset2 = len(dist[list(dist.keys())[0]].keys())
    # get subjects list
    subs_dataset1 = list(sorted(dist.keys()))
    subs_dataset2 = list(sorted(dist[list(dist.keys())[0]].keys()))
    # generate look up table
    gen_lookup_table(subs_dataset1, subs_dataset2, output_path, tp, name1,
                     name2)
    # generate a squared numpy array for storing distance values
    # note we assume nb_subs_dataset2 > nb_subs_dataset1
    cost_matrix = np.zeros((nb_subs_dataset2, nb_subs_dataset2))
    MAX = int(1e6)  # very large number for filling into cost_matrix
    for i in range(nb_subs_dataset1):
        for j in range(nb_subs_dataset2):
            cost_matrix[i, j] = dist[subs_dataset1[i]][
                subs_dataset2[j]]['dist']
    # for the unfilled elements, fill with MAX
    cost_matrix[nb_subs_dataset1:, :].fill(MAX)
    # save cost_matrix
    save_path = os.path.join(output_path, 'cost_matrix_' + str(tp) + '.npy')
    np.save(save_path, cost_matrix)


def gen_lookup_table(subs_dataset1, subs_dataset2, output_path, tp, name1,
                     name2):
    """
    Generate a look up table for searching index in squared cost matrix

    Args:
        subs_dataset1 (list): List of subject RIDs for dataset1
        subs_dataset2 (list): List of subject RIDs for dataset1
        output_path (str): Path for saving output
        tp (int): Number of time points to begin matching
        name1 (str): Name for dataset1
        name2 (str): Name for dataset2
    """
    # columns
    column_names = [name1 + '/' + name2] + \
                   ['j=' + str(j) for j in range(len(subs_dataset2))]
    look_up_lists = list()
    for i, sub_dataset1 in enumerate(subs_dataset1):
        row = ['i=' + str(i)]
        for j, sub_dataset2 in enumerate(subs_dataset2):
            index = str(sub_dataset1) + ',' + str(sub_dataset2)
            row.append(index)
        look_up_lists.append(row)
    # save as a csv file
    save_name = 'lookup_table_' + str(tp) + '.csv'
    list2csv(look_up_lists, column_names, os.path.join(output_path, save_name))


def hangarian(args,
              step3_dist_path,
              output_path,
              threshold,
              name1='MACC',
              name2='ADNI'):
    """
    Hangarian algorithm: https://en.wikipedia.org/wiki/Hungarian_algorithm
    Python implementation: https://software.clapper.org/munkres/index.html
    Args:
        args (tuple): Parameters
        step3_dist_path (str): Path for dist path
        output_path (str): Path for saving output
        threshold (int): Number of time points to begin matching
        name1 (str, optional): Name for dataset1. Defaults to 'MACC'.
        name2 (str, optional): Name for dataset2. Defaults to 'ADNI'.
    """
    # Hungarian calculator
    m = Munkres()
    # Lists for saving matched subjects
    matched_subs_dataset1 = []
    matched_subs_dataset2 = []
    # Run Hungarian matching algorithm for different tps
    for tp in range(threshold, 0, -1):
        cost_matrix = np.load(
            os.path.join(output_path, 'cost_matrix_' + str(tp) + '.npy'))
        # read looking up table
        lookup_df = pd.read_csv(
            os.path.join(output_path, 'lookup_table_' + str(tp) + '.csv'))
        # load cross-dataset distnce pkl
        dist_pkl_name = str(args.matching_pair) + '_' + \
            'AGE' + str(args.age_penalty) + '_' + \
            'SEX' + str(args.sex_penalty) + '_' + \
            'DX' + str(args.dx_penalty) + '_' + \
            'MMSE' + str(args.mmse_penalty) + '_' + \
            'TP' + str(tp) + '.pkl'
        dist = load_pkl(os.path.join(step3_dist_path, dist_pkl_name))
        if len(matched_subs_dataset1) != 0:
            lookup, cost = rm_matched_subs(dist, cost_matrix, lookup_df,
                                           matched_subs_dataset1,
                                           matched_subs_dataset2)
        else:
            cost = cost_matrix
            lookup = lookup_df
        # calculate Hangarian match
        index = m.compute(cost)
        true_index = index[:lookup.shape[0]]
        # find and save matched subjects
        matched_subs_dataset1, matched_subs_dataset2 = \
            find_matched_subs(
                true_index, lookup, dist, tp,
                matched_subs_dataset1, matched_subs_dataset2,
                name1, name2, output_path)


def rm_matched_subs(dist, cost_matrix, lookup_df, matched_subs_dataset1,
                    matched_subs_dataset2):
    """
    Remove matched subjects and re-generate lookup table and dist matrix

    Args:
        dist (dict): Dictionary for distance between dataset1 & dataset2
        cost_matrix (ndarray): Numpy array for global cost
        lookup_df (class DataFrame): Lookup Table
        matched_subs_dataset1 (list): List of matched subject RIDs for dataset1
        matched_subs_dataset2 (list): List of matched subject RIDs for dataset2
    """
    # first deep copy lookup_df and cost_matrix
    _lookup_df = deepcopy(lookup_df)
    _cost_matrix = deepcopy(cost_matrix)
    subs_dataset1 = list(sorted(dist.keys()))
    subs_dataset2 = list(sorted(dist[list(dist.keys())[0]].keys()))
    # get row index for matched pair in lookup_df
    row_index = []
    for sub_dataset1 in matched_subs_dataset1:
        row_index.append(subs_dataset1.index(sub_dataset1))
    # get col index and col names
    col_index = []
    col_names = []
    for sub_dataset2 in matched_subs_dataset2:
        col_index.append(subs_dataset2.index(sub_dataset2))
        col_names.append('j=' + str(subs_dataset2.index(sub_dataset2)))
    # delete rows in lookup_df
    _lookup_df.drop(row_index, axis=0, inplace=True)
    _lookup_df.reset_index(drop=True)
    # delete columns in df
    _lookup_df.drop(col_names, axis=1, inplace=True)
    _lookup_df.reset_index(drop=True)
    # delete rows in cost_matrix
    cost = np.delete(_cost_matrix, row_index, axis=0)
    # delete cols in cost matrix
    cost = np.delete(cost, col_index, axis=1)

    # double check whether the delete operation is correct
    assert cost.shape[0] + len(matched_subs_dataset1) == cost_matrix.shape[0]
    assert cost.shape[1] + len(matched_subs_dataset2) == cost_matrix.shape[1]

    assert _lookup_df.shape[0] + len(
        matched_subs_dataset1) == lookup_df.shape[0]
    assert _lookup_df.shape[1] + len(
        matched_subs_dataset2) == lookup_df.shape[1]

    return _lookup_df, cost


def find_matched_subs(true_index, lookup, dist, tp, matched_subs_dataset1,
                      matched_subs_dataset2, name1, name2, output_path):
    """
    Find matched subject pairs

    Args:
        true_index (list): Index for valid matchign pairs
        lookup (class DataFrame): Lookup Table
        dist (dict): Dictionary for distance between dataset1 & dataset2
        tp (int): Matching time points
        matched_subs_dataset1 (list): List of matched subject RIDs for dataset1
        matched_subs_dataset2 (list): List of matched subject RIDs for dataset2
        name1 (str): Name for dataset1
        name2 (str): Name for dataset2
        output_path (str): Path for saving output
    """
    matched_rows = []
    for idx in true_index:
        matched_row = []
        row_index = idx[0]
        col_index = idx[1] + 1
        strings = lookup.iloc[row_index, col_index].split(',', 2)
        rid_dataset1 = strings[0]
        rid_dataset2 = int(strings[1])
        matched_subs_dataset1.append(rid_dataset1)
        matched_subs_dataset2.append(rid_dataset2)
        dates_dataset1 = list(dist[rid_dataset1][rid_dataset2][name1])
        dates_dataset2 = list(dist[rid_dataset1][rid_dataset2][name2])
        matched_row.append(rid_dataset1)
        matched_row += dates_dataset1
        matched_row.append(rid_dataset2)
        matched_row += dates_dataset2
        matched_rows.append(matched_row)
    # create a column to save
    columns = []
    columns.append(name1 + '_RID')
    for t in range(1, tp + 1):
        tp_name = name1 + '_' + str(t)
        columns.append(tp_name)
    columns.append(name2 + '_RID')
    for t in range(1, tp + 1):
        tp_name = name2 + '_' + str(t)
        columns.append(tp_name)

    save_path = os.path.join(output_path, 'matched_' + str(tp) + '.csv')
    list2csv(matched_rows, columns, save_path)

    return matched_subs_dataset1, matched_subs_dataset2


def merge_matched_tables(output_path, threshold, name1='MACC'):
    """
    Merging matched tables

    Args:
        output_path (str): Path for saving putput
        threshold (int): Threshold to begin matching
        name1 (str, optional): Name for MACC/AIBL. Defaults to 'MACC'.
    """
    matched_df = pd.DataFrame({name1 + '_RID': []})
    for tp in range(threshold, 0, -1):
        matched_csv_path = os.path.join(output_path,
                                        'matched_' + str(tp) + '.csv')
        # load matched csv
        df = pd.read_csv(matched_csv_path)
        matched_df = matched_df.merge(df, how='outer')
    # save to a csv file
    save_path = os.path.join(output_path, 'matched.csv')
    matched_df.to_csv(save_path, index=False, sep=',')
