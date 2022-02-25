#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os


class config:
    NAME = 'ukbb'

    # the files downloaded from original UK Biobank
    REPL_DIR = os.path.join(os.environ['CBIG_REPDATA_DIR'],
                            'stable_projects/predict_phenotypes/He2022_MM')
    CSV_RAW = os.path.join(REPL_DIR, 'df_data_raw.csv')
    DIR_PFC = os.path.join(REPL_DIR, '25753_2_0')

    # config for data processing code
    ROOT_DIR = os.path.join(os.environ['CBIG_CODE_DIR'],
                            'stable_projects/predict_phenotypes/He2022_MM')
    BASE_DIR = os.path.join(ROOT_DIR, 'data_processing')

    # A csv with data-field category information list for UK Biobank
    CSV_UKBB_CAT_INFO = os.path.join(BASE_DIR, 'ukbb_complete_list.csv')

    # config for step 1 coarse filter to phenotypes of UK Biobank dataset
    DIR_1 = os.path.join(BASE_DIR, 'step1_coarse_filter')
    DIR_1_OUTPUT = os.path.join(DIR_1, 'output')

    CSV_BASE = os.path.join(DIR_1_OUTPUT, NAME + '_base.csv')
    CSV_COARSE_FILTER = os.path.join(DIR_1_OUTPUT,
                                     NAME + '_coarse_filter_phe.csv')
    CSV_1000_FOR_PHE_SELECT = os.path.join(DIR_1_OUTPUT,
                                           NAME + '_1000_phe.csv')
    CSV_REMAIN = os.path.join(DIR_1_OUTPUT, NAME + '_remove_1k.csv')
    NPZ_PFC_1000 = os.path.join(DIR_1_OUTPUT, NAME + '_1000_pfc.npz')
    MAT_PFC_1000 = os.path.join(DIR_1_OUTPUT, NAME + '_1000_pfc.mat')
    NPZ_PFC_1000_FLAT = os.path.join(DIR_1_OUTPUT, NAME + '_1000_pfc_flat.npz')
    MAT_PFC_1000_FLAT = os.path.join(DIR_1_OUTPUT, NAME + '_1000_pfc_flat.mat')
    SUBJ_LIST_1000 = os.path.join(DIR_1_OUTPUT, NAME + '_1000_subj_list.npz')
    SUBJ_LIST_1000_TXT = os.path.join(DIR_1_OUTPUT,
                                      NAME + '_1000_subj_list.txt')
    PHE_LIST_COARSE_FILTER = os.path.join(DIR_1_OUTPUT,
                                          NAME + '_coarse_filter_phe_list.txt')

    # config for step 2 KRR to filter phenotypes with relative good prediction
    DIR_2 = os.path.join(BASE_DIR, 'step2_krr_filter')
    DIR_2_OUTPUT = os.path.join(DIR_2, 'output')

    NPZ_PHE_SELECT_RES = os.path.join(DIR_2_OUTPUT, 'phe_filter_result.npz')

    # config for step 3 PCA to part of phenotypes to remove redundant
    DIR_3 = os.path.join(BASE_DIR, 'step3_post_krr_pca')
    DIR_3_OUTPUT = os.path.join(DIR_3, 'output')

    CSV_1000_2ND_FILTER = os.path.join(DIR_3_OUTPUT,
                                       NAME + '_1000_2nd_filter.csv')
    CSV_MAIN_2ND_FILTER = os.path.join(DIR_3_OUTPUT,
                                       NAME + '_main_2nd_filter.csv')
    CSV_TRAIN_2ND_FILTER = os.path.join(DIR_3_OUTPUT,
                                        NAME + '_train_2nd_filter.csv')
    CSV_TEST_2ND_FILTER = os.path.join(DIR_3_OUTPUT,
                                       NAME + '_test_2nd_filter.csv')
    DICT_SAVE_TRAIN = {
        'npz': os.path.join(DIR_3_OUTPUT, NAME + '_train_pfc.npz'),
        'mat': os.path.join(DIR_3_OUTPUT, NAME + '_train_pfc.mat'),
        'npz_flat': os.path.join(DIR_3_OUTPUT, NAME + '_train_pfc_flat.npz'),
        'mat_flat': os.path.join(DIR_3_OUTPUT, NAME + '_train_pfc_flat.mat'),
        'subj_list': os.path.join(DIR_3_OUTPUT, NAME + '_train_subj_list.npz'),
        'subj_list_txt': os.path.join(DIR_3_OUTPUT,
                                      NAME + '_train_subj_list.txt'),
    }
    DICT_SAVE_TEST = {
        'npz': os.path.join(DIR_3_OUTPUT, NAME + '_test_pfc.npz'),
        'mat': os.path.join(DIR_3_OUTPUT, NAME + '_test_pfc.mat'),
        'npz_flat': os.path.join(DIR_3_OUTPUT, NAME + '_test_pfc_flat.npz'),
        'mat_flat': os.path.join(DIR_3_OUTPUT, NAME + '_test_pfc_flat.mat'),
        'subj_list': os.path.join(DIR_3_OUTPUT, NAME + '_test_subj_list.npz'),
        'subj_list_txt': os.path.join(DIR_3_OUTPUT,
                                      NAME + '_test_subj_list.txt'),
    }
    CSV_1000_PCA_PHE = os.path.join(DIR_3_OUTPUT, NAME + '_1000_pca_phe.csv')
    PHE_LIST_PCA = os.path.join(DIR_3_OUTPUT, NAME + '_pca_phe_list.txt')

    DIR_2_OUTPUT_PCA = os.path.join(DIR_2, 'output_pca')
    NPZ_PHE_PCA_SELECT_RES = os.path.join(DIR_2_OUTPUT_PCA,
                                          'phe_pca_select_result.npz')

    # config for step 4 prepare all the data needed for experiment 1
    DIR_4 = os.path.join(BASE_DIR, 'step4_experiment_1')
    DIR_4_OUTPUT = os.path.join(DIR_4, 'output')

    CSV_TRAIN_3RD_FILTER = os.path.join(DIR_4_OUTPUT,
                                        NAME + '_train_3rd_filter.csv')
    CSV_TEST_3RD_FILTER = os.path.join(DIR_4_OUTPUT,
                                       NAME + '_test_3rd_filter.csv')
    CSV_TRAIN_FINAL = os.path.join(DIR_4_OUTPUT, NAME + '_train_final.csv')
    CSV_TEST_FINAL = os.path.join(DIR_4_OUTPUT, NAME + '_test_final.csv')
    PHE_LIST_TRAIN_FINAL = os.path.join(DIR_4_OUTPUT,
                                        NAME + '_train_final_phe_list.txt')
    PHE_LIST_TEST_FINAL = os.path.join(DIR_4_OUTPUT,
                                       NAME + '_test_final_phe_list.txt')
    DICT_TRAIN_TEST_COMBINE = {
        'npz':
        os.path.join(DIR_4_OUTPUT, NAME + '_train_test_pfc.npz'),
        'mat':
        os.path.join(DIR_4_OUTPUT, NAME + '_train_test_pfc.mat'),
        'npz_flat':
        os.path.join(DIR_4_OUTPUT, NAME + '_train_test_pfc_flat.npz'),
        'mat_flat':
        os.path.join(DIR_4_OUTPUT, NAME + '_train_test_pfc_flat.mat'),
        'subj_list':
        os.path.join(DIR_4_OUTPUT, NAME + '_train_test_subj_list.npz'),
        'subj_list_txt':
        os.path.join(DIR_4_OUTPUT, NAME + '_train_test_subj_list.txt'),
        'csv':
        os.path.join(DIR_4_OUTPUT, NAME + '_train_test_final.csv')
    }

    NPZ_FNN_INPUT = os.path.join(DIR_4_OUTPUT, NAME + '_dnn_input_test.npz')

    DIR_FINAL_OUT = os.path.join(BASE_DIR, 'data_processed')

    # config for step 5 HCP data
    DIR_HCP_BASE = '/mnt/isilon/CSC1/Yeolab/Data/HCP/S1200/individuals/'
    DIR_5 = os.path.join(BASE_DIR, 'step5_hcp_data')
    DIR_5_OUTPUT = os.path.join(DIR_5, 'output')
    SUBJ_LIST_HCP_TC = os.path.join(DIR_5_OUTPUT,
                                    'surf_file_list_S1200_1094_210110.txt')
    DIR_5_FC = {
        '419': os.path.join(DIR_5_OUTPUT, 'FC_419'),
    }
    HCP_SUBJ_LIST = os.path.join(DIR_5_OUTPUT, 'HCP_diff_roi_subj_list.txt')
    HCP_PHE_LIST = os.path.join(DIR_5_OUTPUT,
                                'HCP_diff_roi_final_phe_list.txt')

    # config for step 6 prepare all the data needed for experiment 2
    DIR_6 = os.path.join(BASE_DIR, 'step6_experiment_2')
    DIR_6_OUTPUT = os.path.join(DIR_6, 'output')

    DIR_DIFF_ROI_UKBB = {
        '419': os.path.join(DIR_6, 'ukbb_fc_419'),
    }

    CSV_DIFF_ROI_UKBB_2ND_FILTER = os.path.join(
        DIR_6_OUTPUT, NAME + '_diff_roi_2nd_filter.csv')

    DICT_DIFF_ROI = {
        'npz':
        os.path.join(DIR_6_OUTPUT, NAME + '_diff_roi_pfc.npz'),
        'mat':
        os.path.join(DIR_6_OUTPUT, NAME + '_diff_roi_pfc.mat'),
        'npz_flat':
        os.path.join(DIR_6_OUTPUT, NAME + '_diff_roi_pfc_flat.npz'),
        'mat_flat':
        os.path.join(DIR_6_OUTPUT, NAME + '_diff_roi_pfc_flat.mat'),
        'subj_list':
        os.path.join(DIR_6_OUTPUT, NAME + '_diff_roi_subj_list.npz'),
        'subj_list_txt':
        os.path.join(DIR_6_OUTPUT, NAME + '_diff_roi_subj_list.txt'),
    }

    CSV_DIFF_ROI_UKBB_3RD_FILTER = os.path.join(
        DIR_6_OUTPUT, NAME + '_diff_roi_3rd_filter.csv')
    CSV_DIFF_ROI_FINAL = os.path.join(DIR_6_OUTPUT,
                                      NAME + '_diff_roi_final.csv')
    CSV_HCP_DIFF_ROI_FINAL = os.path.join(DIR_6_OUTPUT,
                                          'HCP_diff_roi_final.csv')
    HCP_MAT = os.path.join(DIR_6_OUTPUT, 'HCP_diff_roi_pfc.mat')

    PHE_LIST_DIFF_ROI_FINAL = os.path.join(
        DIR_6_OUTPUT, NAME + '_diff_roi_final_phe_list.txt')

    NPZ_FNN_INPUT_ACROSS_DATASET = os.path.join(
        DIR_6_OUTPUT, NAME + '_dnn_input_cross_dataset.npz')

    # general configs
    RAMDOM_SEED = 42
    KS = [10, 20, 50, 100, 200]

    DICT_PHE_NAME = {
        '12144-2.999': 'Anthropometry Comp.1',
        '12144-2.1000': 'Anthropometry Comp.2',
        '12144-2.1001': 'Anthropometry Comp.3',
        '20009-2.999': 'Non-cancer illness Comp.1',
        '20009-2.1002': 'Non-cancer illness Comp.4',
        '20156-0.999': 'Trail making online Comp.1',
        '20156-0.1002': 'Trail making online Comp.4',
        '20159-0.999': 'Symbol digit substitution online Comp.1',
        '20159-0.1004': 'Symbol digit substitution online Comp.6',
        '21671-2.999': 'Process durations Comp.1',
        '21671-2.1000': 'Process durations Comp.2',
        '21671-2.1002': 'Process durations Comp.4',
        '22022-0.999': 'Genotype Sex inference Comp.1',
        '22022-0.1000': 'Genotype Sex inference Comp.2',
        '30620-0.1000': 'Blood assays Comp.2',
        '30620-0.1001': 'Blood assays Comp.3',
        '30620-0.1002': 'Blood assays Comp.4',
        '30620-0.1003': 'Blood assays Comp.5',
        '3143-0.999': 'Bone-densitometry of heel Comp.1',
        '3143-0.1001': 'Bone-densitometry of heel Comp.3',
        '4194-2.1000': 'Blood pressure & eye measures Comp.2',
        '4194-2.1001': 'Blood pressure & eye measures Comp.3',
        '4194-2.1002': 'Blood pressure & eye measures Comp.4',
        '4194-2.1003': 'Blood pressure & eye measures Comp.5',
        '4194-2.1004': 'Blood pressure & eye measures Comp.6',
        '20007-2.999': 'Age cancer diagnosis Comp.1',
        '31-0.0': 'Sex',
        '3659-0.0': 'Year immigrated to UK',
        '21004-2.999': 'Tower rearranging Comp.1',
        '4286-2.999': 'Prospective memory Comp.1',
        '709-2.0': 'Number in household',
        '26414-0.999': 'Multiple Deprivation Comp.1',
        '20133-0.2': 'Pairs matching online',
        '2946-0.999': 'Family history Comp.1',
        '12336-2.999': 'ECG measures Comp.1',
        '12336-2.1000': 'ECG measures Comp.2',
        '12336-2.1001': 'ECG measures Comp.3',
        '12336-2.1004': 'ECG measures Comp.6',
        '22672-2.999': 'Carotid ultrasound Comp.1',
        '22672-2.1003': 'Carotid ultrasound Comp.5',
        '6373-2.999': 'Matrix pattern completion Comp.1',
        '6373-2.1000': 'Matrix pattern completion Comp.2',
        '6373-2.1001': 'Matrix pattern completion Comp.3',
        '20162-2.999': 'Smoke Comp.1',
        '46-2.999': 'Hand grip strength Comp.1',
        '23323-2.999': 'Symbol digit substitution Comp.1',
        '4440-2.0': 'Monthly spirits',
        '20127-0.0': 'Neuroticism score',
        '20075-2.999': 'Location Comp.1',
        '767-2.0': 'Weekly working hour ',
        '20016-2.0': 'Fluid intelligence',
        '4230-2.3': 'Hearing SNR test',
        '1070-2.0': 'Time spent watching TV',
        '864-2.0': '#days/week walked 10+ minutes',
        '22003-0.999': 'Genetic Comp.1',
        '3062-2.999': 'Spirometry Comp.1',
        '1160-2.0': 'Sleep duration',
        '845-2.0': 'Age completed full time education',
        '6348-2.999': 'Trail making Comp.1',
        '398-2.3': 'Pairs matching',
        '1578-2.0': 'Weekly champagne + white wine',
        '4250-2.999': 'Numeric memory Comp.1',
        '1588-2.0': 'Weekly beer plus cider',
        '777-2.0': 'Frequency home to workplace',
        'age-2.0': 'Age',
        '30512-1.999': 'Urine assays Comp.1',
        '1090-2.0': 'Time spent driving',
        '20023-2.999': 'Reaction time Comp.1',
    }
