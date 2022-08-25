#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import logging
import numpy as np
import pandas as pd
from utils.misc import create_folder


def split_MACC_into_bins(args, BINS=[]):
    """
    Split MACC subjects into bins according to subject-level MMSE

    Args:
        args (tuple): Parameters
        BINS (list, optional): BINs to split MACC subjects into.
    """
    # set output path
    output_path = os.path.join(args.checkpoint_path, args.matching_pair,
                               'matching_' + str(args.nb_bins) + 'BINs')
    create_folder(output_path)
    create_folder(os.path.join(output_path, 'MERGED'))
    # read raw data
    raw_dataset1 = pd.read_csv(args.MACC_data_path)
    # calculate subject-level MMSE
    subjects = np.unique(raw_dataset1.RID)
    nb_subjects = len(subjects)
    rows = []
    cols = ['RID', 'sMMSE']
    for sub in subjects:
        sub_mask = (raw_dataset1.RID == sub)
        sub_data = raw_dataset1[sub_mask]
        sub_row = []
        sub_row = list(sub_data[['MMSE']].mean().values)
        sub_row = [sub] + sub_row
        rows.append(sub_row)
    dataset1_subject = pd.DataFrame(data=rows, columns=cols)
    # equally set bins or manually set bins
    if len(BINS) == 0:
        MMSEs = np.arange(31.0)  # MMSE ranges from 0 to 30
        MMSE_bins = np.array_split(MMSEs, args.nb_bins)
    else:
        MMSE_bins = BINS
    _subjects = list(subjects)
    for bin in range(args.nb_bins):
        bin_output_path = os.path.join(output_path, 'BIN_' + str(bin))
        create_folder(bin_output_path)
        MMSE_bin = MMSE_bins[bin]
        subject_df = dataset1_subject[(dataset1_subject.sMMSE >= MMSE_bin[0]) &
                                      (dataset1_subject.sMMSE < MMSE_bin[-1])]
        subs = list(np.unique(subject_df.RID))
        _subjects = list(set(_subjects) - set(subs))
        if bin == args.nb_bins - 1:
            subs = subs + _subjects
            subject_df = subject_df.append(
                dataset1_subject[(dataset1_subject.sMMSE == MMSE_bin[-1])],
                ignore_index=True)
        df = pd.DataFrame(columns=list(raw_dataset1.columns))
        for sub in subs:
            sub_mask = (raw_dataset1.RID == sub)
            sub_data = raw_dataset1[sub_mask]
            df = df.append(sub_data, ignore_index=True)
        # study data distribution
        study_bins_distribution(df, subject_df, nb_subjects, MMSE_bin,
                                args.nb_bins, bin, bin_output_path)
        # save as csv files
        df.to_csv(
            os.path.join(bin_output_path, 'MACC_' + str(bin) + '_bin.csv'),
            sep=',',
            index=False)
    create_folder(
        os.path.join(output_path, 'BIN_' + str(args.nb_bins), 'round_1'))


def study_bins_distribution(df, subject_df, nb_subjects, MMSE_bin, nb_bins,
                            bin, output_path):
    """
    1. Calculate #TPs and #Subjects within this bin
    2. Calculate #TPs and #Subjects for each MMSE score (0, 1, 2.....)
       and write to a log file

    Args:
        df (class Dataframe): MACC dataframe
        subject_df (class DataFrame): Subject dataframe
        nb_subjects (int): Number of subjects
        MMSE_bin (list): List of bins
        nb_bins (int): Number of bins
        bin (int): Bin
        output_path (str): Path for saving output
    """
    # config logging module
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    fh = logging.FileHandler(
        filename=os.path.join(output_path, 'BIN_' + str(bin) + '_info.txt'))
    fh.setLevel(logging.INFO)
    log.addHandler(fh)

    log.info('Total #subjects in MACC is ' + str(nb_subjects))
    log.info('This is ' + str(bin) + 'th bin out of ' + str(nb_bins) + ' BINs')
    log.info('Total #subjects within this bin:' +
             str(np.unique(df.RID).shape[0]))
    log.info('Total #TPs within this bin:' + str(df.shape[0]))
    # for each MMSE score, get how many subjects and TPs
    log.info('Detalied Subject-level info within this bin..')
    for MMSE_score in MMSE_bin:
        MMSE_df = subject_df[(subject_df.sMMSE > MMSE_score - 1)
                             & (subject_df.sMMSE <= MMSE_score)]
        nb_subs_MMSE_df = np.unique(MMSE_df.RID).shape[0]
        log.info('...MMSE <=' + str(MMSE_score) + '; #Subjects=' +
                 str(nb_subs_MMSE_df) + '; #TPs' + str(MMSE_df.shape[0]))
    log.removeHandler(fh)
    del log, fh
