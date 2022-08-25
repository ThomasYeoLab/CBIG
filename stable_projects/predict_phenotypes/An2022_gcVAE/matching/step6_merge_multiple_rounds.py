#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import pandas as pd
from utils.misc import create_folder


def merge_multiple_rounds(args,
                          bin,
                          round1_threshold,
                          round2_threshold,
                          name1='MACC',
                          name2='ADNI'):
    """
    Merge data (picked_adni.csv & picked_macc.csv) from multiple rounds

    Args:
        args (tuple): Parameters
        bin (int): Bin
        round1_threshold (int): Threshold for round 1 matching
        round2_threshold (int): Threshold for round 2 matching
        name1 (str, optional): Name for dataset1. Defaults to 'MACC'.
        name2 (str, optional): Name for dataset2. Defaults to 'ADNI'.
    """
    # read round1 data
    round1_step5_output_path = os.path.join(
        args.checkpoint_path, args.matching_pair,
        'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin), 'round_1',
        'step5_picking',
        'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
        + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty),
        'Threshold_' + str(round1_threshold))
    round1_macthed_dataset1 = pd.read_csv(
        os.path.join(round1_step5_output_path, 'picked_' + name1 + '.csv'))
    round1_macthed_dataset2 = pd.read_csv(
        os.path.join(round1_step5_output_path, 'picked_' + name2 + '.csv'))

    # read round 2 data
    threshold_name = 'Threshold_' + str(round1_threshold) + '-' + str(
        round2_threshold)
    round2_step5_output_path = os.path.join(
        args.checkpoint_path, args.matching_pair,
        'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin), 'round_2',
        'step5_picking',
        'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
        + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty),
        'Threshold_' + str(round2_threshold))
    round2_macthed_dataset1 = pd.read_csv(
        os.path.join(round2_step5_output_path, 'picked_' + name1 + '.csv'))
    round2_macthed_dataset2 = pd.read_csv(
        os.path.join(round2_step5_output_path, 'picked_' + name2 + '.csv'))

    merged_dataset1 = round1_macthed_dataset1.append(
        round2_macthed_dataset1, sort=False)
    merged_dataset2 = round1_macthed_dataset2.append(
        round2_macthed_dataset2, sort=False)

    output_path = os.path.join(
        args.checkpoint_path, args.matching_pair,
        'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin), 'merged',
        'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
        + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty),
        threshold_name)
    create_folder(output_path)
    # save merged data
    merged_dataset1.to_csv(
        os.path.join(output_path, 'picked_' + name1 + '.csv'),
        sep=',',
        index=False)
    merged_dataset2.to_csv(
        os.path.join(output_path, 'picked_' + name2 + '.csv'),
        sep=',',
        index=False)
