#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import pandas as pd
import numpy as np
from copy import deepcopy
from utils.misc import load_pkl, save_pkl
from matching.step5_pick_well_matched_pairs\
    import gen_matched_subjects_cost_table, rm_high_cost_pairs


def control_within_bin_matched_pairs(args, round1_threshold, round2_threshold,
                                     bin):
    """
    Control #matched pairs within each bin

    Args:
        args (tuple): Parameters
        round1_threshold (int): Threshold for round 1 matching
        round2_threshold (int): Threshold for round 2 matching
        bin (int): Bin
    """
    curr_bin_path = os.path.join(args.checkpoint_path, args.matching_pair,
                                 'matching_' + str(args.nb_bins) + 'BINs',
                                 'BIN_' + str(bin))
    curr_merged_path = os.path.join(
        curr_bin_path, 'merged',
        'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
        + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty),
        'Threshold_' + str(round1_threshold) + '-' + str(round2_threshold))
    # read raw dataset
    raw_dataset1 = pd.read_csv(
        os.path.join(args.checkpoint_path, args.matching_pair,
                     'matching_' + str(args.nb_bins) + 'BINs',
                     'BIN_' + str(bin), 'MACC_' + str(bin) + '_bin.csv'))
    nb_bin_subjects = len(np.unique(raw_dataset1.RID))
    raw_dataset2 = pd.read_csv(args.ADNI_data_path)
    # read picked_matched
    picked_MACC = pd.read_csv(
        os.path.join(curr_merged_path, 'picked_MACC.csv'))
    picked_ADNI = pd.read_csv(
        os.path.join(curr_merged_path, 'picked_ADNI.csv'))

    assert picked_ADNI.shape[0] == picked_MACC.shape[0], 'Wrong matched case'

    if picked_MACC.shape[0] <= args.match_ratio * nb_bin_subjects:
        # no need to control
        picked_MACC.to_csv(
            os.path.join(curr_merged_path, 'controled_MACC.csv'),
            sep=',',
            index=False)
        picked_ADNI.to_csv(
            os.path.join(curr_merged_path, 'controled_ADNI.csv'),
            sep=',',
            index=False)
    else:
        # we need to drop some subjects according to matching cost
        # note that do not let DX==0 have too much for 1tp subject
        matched_cost = gen_matched_subjects_cost_table(
            args, raw_dataset1, raw_dataset2, picked_MACC, picked_ADNI,
            curr_merged_path)
        rm_macc_subs = matched_cost.iloc[int(args.match_ratio *
                                             nb_bin_subjects):, 0].values
        rm_adni_subs = matched_cost.iloc[int(args.match_ratio *
                                             nb_bin_subjects):, 1].values
        # remove rm_macc_subs from picked_MACC.csv
        rm_high_cost_pairs(
            rm_macc_subs, picked_MACC,
            os.path.join(curr_merged_path, 'controled_MACC.csv'))
        # remove rm_adni_subs from picked_MACC.csv
        rm_high_cost_pairs(
            rm_adni_subs, picked_ADNI,
            os.path.join(curr_merged_path, 'controled_ADNI.csv'))


def update_next_bin_data(args,
                         bin,
                         round1_threshold,
                         round2_threshold,
                         name1='MACC',
                         name2='ADNI'):
    """
    Using matched ADNI subjects to update next bin ADNI.pkl

    Args:
        args (tuple): Parameters
        bin (int): Bin
        round1_threshold (int): Threshold for round 1 matching
        round2_threshold (int): Threshold for round 2 matching
        name1 (str, optional): Name for dataset1. Defaults to 'MACC'.
        name2 (str, optional): Name for dataset2. Defaults to 'ADNI'.
    """
    curr_bin_path = os.path.join(args.checkpoint_path, args.matching_pair,
                                 'matching_' + str(args.nb_bins) + 'BINs',
                                 'BIN_' + str(bin))
    next_bin_path = os.path.join(args.checkpoint_path, args.matching_pair,
                                 'matching_' + str(args.nb_bins) + 'BINs',
                                 'BIN_' + str(bin + 1))
    # read last round ADNI data in current bin
    curr_last_round_path = os.path.join(curr_bin_path,
                                        'round_' + str(args.nb_rounds + 1))
    ADNI_pkl_name = args.name_2 + '_' + 'AGE' + str(args.age_penalty) + '_' + \
        'SEX' + str(args.sex_penalty) + '_' + \
        'DX' + str(args.dx_penalty) + '_' +\
        'MMSE' + str(args.mmse_penalty) + '.pkl'
    ADNI_dist = load_pkl(os.path.join(curr_last_round_path, ADNI_pkl_name))
    curr_merged_path = os.path.join(
        curr_bin_path, 'merged',
        'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
        + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty),
        'Threshold_' + str(round1_threshold) + '-' + str(round2_threshold))
    picked_ADNI = pd.read_csv(
        os.path.join(curr_merged_path, 'controled_ADNI.csv'))
    matched_subs = np.unique(picked_ADNI.RID)
    _ADNI_dist = deepcopy(ADNI_dist)
    for sub in _ADNI_dist.keys():
        if sub in matched_subs:
            del ADNI_dist[sub]
    # save to next bin round1
    output_path = os.path.join(next_bin_path, 'round_1', ADNI_pkl_name)
    save_pkl(ADNI_dist, output_path)
