#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import shutil
import warnings
import argparse
from config import global_config
from utils.misc import\
    create_folder, myException, load_pkl, clean_folder
# import functions for matching
from matching.matching_misc\
    import get_max_threshold_round, get_max_threshold
from matching.step1_split_subjects_into_bins\
    import split_MACC_into_bins
from matching.step2_get_within_dataset_combs\
    import get_within_dataset_combs
from matching.step3_get_cross_dataset_dist\
    import get_cross_dataset_dist
from matching.step4_run_Hangarian_matching\
    import hangarian_matching_wrapper
from matching.step5_pick_well_matched_pairs\
    import pick_well_matched_pairs, split_matched_csv, update_next_round_data
from matching.step6_merge_multiple_rounds\
    import merge_multiple_rounds
from matching.step7_control_matched_quantity_within_bins\
    import control_within_bin_matched_pairs, update_next_bin_data
from matching.step8_merge_multiple_bins\
    import merge_multi_bins, extract_matched_unmatched_data


def matching_args_parser():
    """
    Arugment parser for matching ADNI-MACC (ADNI-AIBL)

    """
    parser = argparse.ArgumentParser(prog='MatchingArgs')
    parser.add_argument(
        '--ADNI_data_path',
        type=str,
        default=os.path.join(global_config.raw_data_path, 'ADNI.csv'))
    parser.add_argument(
        '--MACC_data_path',
        type=str,
        default=os.path.join(global_config.raw_data_path, 'MACC.csv'))
    parser.add_argument(
        '--columns_path', type=str, default=global_config.columns_path)
    parser.add_argument(
        '--checkpoint_path',
        type=str,
        default=os.path.join(global_config.checkpoints_path, 'matching'))
    parser.add_argument('--matching_pair', '-p', type=str, default='ADNI-MACC')
    parser.add_argument('--name_1', type=str, default='MACC')
    parser.add_argument('--name_2', type=str, default='ADNI')

    parser.add_argument('--nb_bins', type=int, default=3)
    parser.add_argument('--nb_rounds', type=int, default=2)
    parser.add_argument('--match_ratio', type=float, default=0.5)
    parser.add_argument('--round', type=int, default=1)
    parser.add_argument('--threshold', type=int, default=4)
    parser.add_argument('--NANpenalty', type=float, default=10.0)
    parser.add_argument('--sex_penalty', type=float, default=10.0)
    parser.add_argument('--age_penalty', type=float, default=5.0)
    parser.add_argument('--dx_penalty', type=float, default=10.0)
    parser.add_argument('--mmse_penalty', type=float, default=10.0)
    parser.add_argument('--mmse_eps', type=float, default=2.0)

    args, _ = parser.parse_known_args()

    return args


def matching_wrapper(args):
    """
    Top level wrapper for matching subjects between ADNI-MACC (ADNI-AIBL)

    Args:
        args (tuple): Parameters
    """
    # check for ADNI-MACC or ADNI-AIBL
    if args.matching_pair == 'ADNI-MACC':
        args.MACC_data_path = os.path.join(global_config.raw_data_path,
                                           'MACC.csv')
        args.name1 = 'MACC'
    elif args.matching_pair == 'ADNI-AIBL':
        args.MACC_data_path = os.path.join(global_config.raw_data_path,
                                           'AIBL.csv')
        args.name1 = 'AIBL'
    else:
        raise myException(
            'args.matching_pair can only be ADNI-AIBL or ADNI-MACC')
    # step 1
    matching_step1_wrapper(args)
    # step 2
    matching_step2_wrapper(args)
    # step 3, 4, 5, 6, 7
    for bin in range(args.nb_bins):
        matching_step345_wrapper(args, bin)
        matching_step6_wrapper(args, bin)
        matching_step7_wrapper(args, bin)
    # step 8
    matching_step8_wrapper(args)
    # copy matched and unmatched data to "An2022_gcVAE/data/splits" folder
    src_csvs = [
        'ADNI_matched.csv', 'ADNI_unmatched.csv', 'MACC_matched.csv',
        'MACC_unmatched.csv'
    ]
    src_dir = os.path.join(args.checkpoint_path, args.matching_pair,
                           'matching_3BINs', 'MERGED',
                           'AGE5.0_SEX10.0_DX10.0_MMSE10.0')
    dst_csvs = [
        'ADNI_matched.csv', 'ADNI_unmatched.csv', args.name1 + '_matched.csv',
        args.name1 + '_unmatched.csv'
    ]
    dst_dir = os.path.join(global_config.data_path, 'splits',
                           args.matching_pair)
    create_folder(dst_dir)
    for src_csv, dst_csv in zip(src_csvs, dst_csvs):
        shutil.copyfile(
            os.path.join(src_dir, src_csv), os.path.join(dst_dir, dst_csv))
    # clean checkpoints folder
    clean_folder(os.path.join(args.checkpoint_path, args.matching_pair))


def matching_step1_wrapper(args):
    """
    Wrapper function for step1_split_subjects_into_bins

    Args:
        args (tuple): Parameters
    """
    split_MACC_into_bins(args)


def matching_step2_wrapper(args):
    """
    Wrapper function for step2_get_within_dataset_combs

    Args:
        args (tuple): Parameters
    """
    for bin in range(args.nb_bins):
        if bin == 0:
            # why we only compute withincombs for ADNI only at 0th bin
            # becasue we will update the combs for following bins
            ADNI_output_path = os.path.join(
                args.checkpoint_path, args.matching_pair,
                'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
                'round_' + str(args.round), 'ADNI.pkl')
            create_folder(
                os.path.join(args.checkpoint_path, args.matching_pair,
                             'matching_' + str(args.nb_bins) + 'BINs',
                             'BIN_' + str(bin), 'round_' + str(args.round)))
            get_within_dataset_combs(args.ADNI_data_path, ADNI_output_path)
        # compute for MACC/AIBL
        MACC_input_path = os.path.join(
            args.checkpoint_path, args.matching_pair,
            'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
            'MACC_' + str(bin) + '_bin.csv')
        MACC_output_path = os.path.join(
            args.checkpoint_path, args.matching_pair,
            'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
            'round_' + str(args.round), 'MACC.pkl')
        create_folder(
            os.path.join(args.checkpoint_path, args.matching_pair,
                         'matching_' + str(args.nb_bins) + 'BINs',
                         'BIN_' + str(bin), 'round_' + str(args.round)))
        get_within_dataset_combs(MACC_input_path, MACC_output_path)


def matching_step345_wrapper(args, bin):
    """
    Wrapper function for steps:
        step3_get_cross_dataset_dist
        step4_run_Harngarian_matching
        step5_pick_well_matched_pairs

    Args:
        args (tuple): Parameters
        bin (int): Bin
    """
    eps_list = [10, 3, 2]
    args.mmse_eps = eps_list[bin]
    for round in range(1, args.nb_rounds + 1):
        args.round = round
        # step3_get_cross_dataset_dist
        save_name1 = args.name_1 + '_' + 'AGE' + \
            str(args.age_penalty) + '_' + \
            'SEX' + str(args.sex_penalty) + '_' + \
            'DX' + str(args.dx_penalty) + '_' + \
            'MMSE' + str(args.mmse_penalty) + '.pkl'
        save_name2 = args.name_2 + '_' + 'AGE' + \
            str(args.age_penalty) + '_' + \
            'SEX' + str(args.sex_penalty) + '_' + \
            'DX' + str(args.dx_penalty) + '_' +\
            'MMSE' + str(args.mmse_penalty) + '.pkl'
        step3_output_path = os.path.join(
            args.checkpoint_path, args.matching_pair,
            'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin))
        # load data
        if args.round == 1 and bin == 0:
            adni_data = load_pkl(
                os.path.join(step3_output_path, 'round_' + str(round),
                             'ADNI.pkl'))
            macc_data = load_pkl(
                os.path.join(step3_output_path, 'round_' + str(round),
                             'MACC.pkl'))
        elif args.round == 1 and bin > 0:
            adni_data = load_pkl(
                os.path.join(step3_output_path, 'round_' + str(round),
                             save_name2))
            macc_data = load_pkl(
                os.path.join(step3_output_path, 'round_' + str(round),
                             'MACC.pkl'))
        else:
            adni_data = load_pkl(
                os.path.join(step3_output_path, 'round_' + str(round),
                             save_name2))
            macc_data = load_pkl(
                os.path.join(step3_output_path, 'round_' + str(round),
                             save_name1))
        # get maximun threshold for this round
        max_threshold_round = get_max_threshold_round(macc_data)
        for threshold in range(1, max_threshold_round + 1):
            args.threshold = threshold
            get_cross_dataset_dist(args, bin, macc_data, adni_data)
        # step4_run_Harngarian_matching
        for threshold in range(1, max_threshold_round + 1):
            args.threshold = threshold
            step4_output_path = os.path.join(
                args.checkpoint_path, args.matching_pair,
                'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
                'round_' + str(args.round), 'step4_matching', 'AGE' + str(
                    args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
                + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty))
            create_folder(step3_output_path)
            hangarian_matching_wrapper(args, bin, step4_output_path)
        # step5_pick_well_matched_pairs
        for threshold in range(1, max_threshold_round + 1):
            args.threshold = threshold
            step4_output_path = os.path.join(
                args.checkpoint_path, args.matching_pair,
                'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
                'round_' + str(args.round), 'step4_matching', 'AGE' + str(
                    args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
                + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty))
            step5_output_path = os.path.join(
                args.checkpoint_path, args.matching_pair,
                'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin),
                'round_' + str(args.round), 'step5_picking', 'AGE' + str(
                    args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
                + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty))
            create_folder(step4_output_path)
            pick_well_matched_pairs(args, step4_output_path, step5_output_path)
            matched_subs_dataset1, matched_subs_dataset2 = split_matched_csv(
                args, bin, step5_output_path)
            update_next_round_data(args, bin, matched_subs_dataset1,
                                   matched_subs_dataset2)


def matching_step6_wrapper(args, bin):
    """
    Wrapper function for step6_merge_multiple_rounds

    Args:
        args (tuple): Paramters
        bin (int): Bin
    """
    # note that we need to get the round1_threshold
    round1_threshold = get_max_threshold(
        os.path.join(args.checkpoint_path, args.matching_pair,
                     'matching_' + str(args.nb_bins) + 'BINs',
                     'BIN_' + str(bin), 'MACC_' + str(bin) + '_bin.csv'))
    save_name1 = args.name_1 + '_' + 'AGE' + str(args.age_penalty) + '_' + \
        'SEX' + str(args.sex_penalty) + '_' + \
        'DX' + str(args.dx_penalty) + '_' +\
        'MMSE' + str(args.mmse_penalty) + '.pkl'
    step2_output_path = os.path.join(args.checkpoint_path, args.matching_pair,
                                     'matching_' + str(args.nb_bins) + 'BINs',
                                     'BIN_' + str(bin))
    macc_data = load_pkl(
        os.path.join(step2_output_path, 'round_' + str(2), save_name1))
    round2_threshold = get_max_threshold_round(macc_data)

    merge_multiple_rounds(args, bin, round1_threshold, round2_threshold)


def matching_step7_wrapper(args, bin):
    """
    Wrapper function for step7_control_matched_quantity_within_bins

    Args:
        args (tuple): Paramters
        bin (int): Bin
    """
    round1_threshold = get_max_threshold(
        os.path.join(args.checkpoint_path, args.matching_pair,
                     'matching_' + str(args.nb_bins) + 'BINs',
                     'BIN_' + str(bin), 'MACC_' + str(bin) + '_bin.csv'))
    save_name1 = args.name_1 + '_' + 'AGE' + str(args.age_penalty) + '_' + \
        'SEX' + str(args.sex_penalty) + '_' + \
        'DX' + str(args.dx_penalty) + '_' +\
        'MMSE' + str(args.mmse_penalty) + '.pkl'
    step2_output_path = os.path.join(args.checkpoint_path, args.matching_pair,
                                     'matching_' + str(args.nb_bins) + 'BINs',
                                     'BIN_' + str(bin))
    macc_data = load_pkl(
        os.path.join(step2_output_path, 'round_' + str(2), save_name1))
    round2_threshold = get_max_threshold_round(macc_data)
    control_within_bin_matched_pairs(args, round1_threshold, round2_threshold,
                                     bin)
    update_next_bin_data(args, bin, round1_threshold, round2_threshold)


def matching_step8_wrapper(args):
    """
    Wrapper function for step8_merge_multiple_bins

    Args:
        args (tuple): Parameters
    """
    # get round1_threshold_lists
    round1_threshold_list = []
    round2_threshold_list = []
    for bin in range(args.nb_bins):
        round1_threshold = get_max_threshold(
            os.path.join(args.checkpoint_path, args.matching_pair,
                         'matching_' + str(args.nb_bins) + 'BINs',
                         'BIN_' + str(bin), 'MACC_' + str(bin) + '_bin.csv'))
        round1_threshold_list.append(round1_threshold)
        save_name1 = args.name_1 + '_' + 'AGE' + \
            str(args.age_penalty) + '_' + \
            'SEX' + str(args.sex_penalty) + '_' + \
            'DX' + str(args.dx_penalty) + '_' +\
            'MMSE' + str(args.mmse_penalty) + '.pkl'
        step2_output_path = os.path.join(
            args.checkpoint_path, args.matching_pair,
            'matching_' + str(args.nb_bins) + 'BINs', 'BIN_' + str(bin))
        macc_data = load_pkl(
            os.path.join(step2_output_path, 'round_' + str(2), save_name1))
        round2_threshold = get_max_threshold_round(macc_data)
        round2_threshold_list.append(round2_threshold)

    merge_multi_bins(args, round1_threshold_list, round2_threshold_list)
    # extract matched and unmatched data
    bins_path = os.path.join(args.checkpoint_path, args.matching_pair,
                             'matching_' + str(args.nb_bins) + 'BINs')
    merged_path = os.path.join(
        bins_path, 'MERGED',
        'AGE' + str(args.age_penalty) + '_SEX' + str(args.sex_penalty) + '_DX'
        + str(args.dx_penalty) + '_MMSE' + str(args.mmse_penalty))
    create_folder(merged_path)
    extract_matched_unmatched_data(args, merged_path)


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    matching_wrapper(matching_args_parser())
