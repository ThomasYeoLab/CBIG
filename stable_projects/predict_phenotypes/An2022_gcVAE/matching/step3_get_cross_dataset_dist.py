#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import numpy as np
from utils.misc import save_pkl, create_folder


def get_cross_dataset_dist(args, bin, dataset1, dataset2):
    """
    Compute Cross datset distance by giving parameters.
    Note we are matching subjects in dataset1 to subjects in dataset2,
    which means we assume that dataset1 is smaller than dataset2

    Args:
        args (tuple): Parameters
        bin (int): Bin
        dataset1 (dict): Dictonary for the combinations of dataset1
        dataset2 (dict): Dictonary for the combinations of dataset1
    """
    name1 = args.name_1
    name2 = args.name_2
    dist_matrix = dict()
    subjects_dataset1 = sorted(dataset1.keys())
    subjects_dataset2 = sorted(dataset2.keys())
    for sub_dataset1 in subjects_dataset1:
        # for all possible subject pairs for given threshold
        subdata_dataset1 = dataset1[sub_dataset1]
        if subdata_dataset1['NumTPs'] >= args.threshold:
            dist_matrix[sub_dataset1] = dict()
            for sub_dataset2 in subjects_dataset2:
                subdata_dataset2 = dataset2[sub_dataset2]
                if subdata_dataset2['NumTPs'] >= args.threshold:
                    # calculate minimal distance between 2 datasets
                    dist_matrix[sub_dataset1][sub_dataset2] =\
                         find_min_dist(args, args.threshold,
                                       subdata_dataset1, subdata_dataset2,
                                       name1, name2)
    # save distance matrix
    save_name = str(args.matching_pair) + '_' + \
        'AGE' + str(args.age_penalty) + '_' + \
        'SEX' + str(args.sex_penalty) + '_' + \
        'DX' + str(args.dx_penalty) + '_' + \
        'MMSE' + str(args.mmse_penalty) + '_' + \
        'TP' + str(args.threshold) + '.pkl'
    save_dir = os.path.join(
        args.checkpoint_path, args.matching_pair,
        'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
        'round_' + str(args.round), 'step3_distance',
        'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
        + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty))
    create_folder(save_dir)
    save_path = os.path.join(save_dir, save_name)
    # save to a .pkl file
    save_pkl(dist_matrix, save_path)


def find_min_dist(args, nb_tps, subdata_dataset1, subdata_dataset2, name1,
                  name2):
    """
    Find minmial distance for 2 subjects cross datasets

    Args:
        args (tuple): Parameters
        nb_tps (int): Number of time points
        subdata_dataset1 (dict): Dictonary of subject from dataset1
        subdata_dataset2 (dict): Dictonary of subject from dataset1
        name1 (str): Name for dataset1
        name2 (str): Name for dataset2
    """
    ret = {}
    # calculate distance vectos for all measures
    age_dist = one_measure_dist(args.age_penalty, args.NANpenalty,
                                subdata_dataset1['AGE'][nb_tps]['Matrix'],
                                subdata_dataset2['AGE'][nb_tps]['Matrix'])
    dx_dist = one_measure_dist(args.dx_penalty, args.NANpenalty,
                               subdata_dataset1['DX'][nb_tps]['Matrix'],
                               subdata_dataset2['DX'][nb_tps]['Matrix'])
    mmse_dist = one_measure_dist(args.mmse_penalty, args.NANpenalty,
                                 subdata_dataset1['MMSE'][nb_tps]['Matrix'],
                                 subdata_dataset2['MMSE'][nb_tps]['Matrix'])
    if subdata_dataset1['SEX'] == subdata_dataset2['SEX']:
        sex_dist_value = 0
    else:
        sex_dist_value = args.sex_penalty
    # find minmum distance
    min_dist = age_dist + dx_dist + mmse_dist
    min_dist_index = np.argmin(min_dist)
    min_dist_value = min_dist[min_dist_index] + sex_dist_value
    ret['dist'] = min_dist_value

    # find the corresponding combinatin in 2 datsets
    subage_dataset2 = subdata_dataset2['AGE'][nb_tps]['Matrix']
    subcomb_dataset1 = \
        subdata_dataset1['AGE'][nb_tps]['Comb'][
            min_dist_index // subage_dataset2.shape[0]]

    subcomb_dataset2 = \
        subdata_dataset2['AGE'][nb_tps]['Comb'][
            min_dist_index % subage_dataset2.shape[0]]

    # get corresponding dates
    dates_dataset1 = []
    dates_dataset2 = []
    for tp in range(nb_tps):
        dates_dataset1.append(subdata_dataset1['TPs'][subcomb_dataset1[tp]][0])
        dates_dataset2.append(subdata_dataset2['TPs'][subcomb_dataset2[tp]][0])

    ret[name1] = dates_dataset1
    ret[name2] = dates_dataset2

    return ret


def one_measure_dist(penalty, NANpenalty, matrix_dataset1, matrix_dataset2):
    """
    Find distance for 2 subjects cross datasets on 1 measure

    Args:
        penalty (int): Penalty for penishing badly matched
        NANpenalty (int): Penalty for penishing missing data
        matrix_dataset1 (ndarray): Numpy array for dataset1
        matrix_dataset2 (ndarry): Numpy array for dataset2
    """
    repeated_matrix_dataset1 = np.repeat(
        matrix_dataset1, matrix_dataset2.shape[0], axis=0)
    repeated_matrix_dataset2 = np.vstack(
        [matrix_dataset2] * matrix_dataset1.shape[0])
    dist = np.abs(repeated_matrix_dataset1 -
                  repeated_matrix_dataset2) * penalty

    # replace Nan value
    dist[np.isnan(dist)] = NANpenalty
    dist = np.sum(dist, axis=1)

    return dist
