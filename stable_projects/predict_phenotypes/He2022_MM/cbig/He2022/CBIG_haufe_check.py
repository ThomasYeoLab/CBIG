#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import time
import random
import argparse
import itertools
import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
from scipy.stats import ttest_rel
from config import config
import matplotlib.pyplot as plt
import seaborn as sns
plt.switch_backend('Agg')


def sum_of_mul(A, B):
    '''sum of multiplication of two array over axis=1

    Args:
        A (ndarray): first array for calculation
        B (ndarray): second array for calculation

    Returns:
        ndarray: sum of multiplication calculated
    '''
    return np.einsum('ij,ij->i', A, B)


def covariance_rowwise(A, B):
    '''compute rowwise covariance

    Args:
        A (ndarray): first array for covariance calculation
        B (ndarray): second array for covariance calculation

    Returns:
        ndarray: covariance calculated
    '''
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(0, keepdims=True)
    B_mB = B - B.mean(0, keepdims=True)

    N = A.shape[0]
    if B_mB.ndim == 1:
        B_mB = np.expand_dims(B_mB, -1)
    a_nsample = A_mA.shape[1]
    b_nsample = B_mB.shape[1]
    rnt = np.zeros((a_nsample, b_nsample))
    comb = np.array(
        list(itertools.product(range(a_nsample), range(b_nsample))))

    n_comb = len(comb)
    chunk = 100000
    if n_comb > chunk:
        start_time = time.time()
        cov = np.empty(n_comb)
        for i in range(chunk, n_comb, chunk):
            cov[i - chunk:i] = sum_of_mul(A_mA[:, comb[i - chunk:i, 0]].T,
                                          B_mB[:, comb[i - chunk:i, 1]].T)

            print(i, time.time() - start_time)
        cov[i:] = sum_of_mul(A_mA[:, comb[i:, 0]].T, B_mB[:, comb[i:, 1]].T)
    else:
        cov = sum_of_mul(A_mA[:, comb[:, 0]].T, B_mB[:, comb[:, 1]].T)
    rnt[comb[:, 0], comb[:, 1]] = cov
    return np.squeeze(rnt) / (N - 1)


def compute_PNF(x, y_pred):
    '''compute predictive network features (PNF) with various shape

    Args:
        x (ndarray): array of FC for PNF calculation
        y_pred (ndarray): array of predicted y for PNF calculation

    Returns:
        ndarray: PNF calculated
    '''
    if len(x.shape) < 3 and len(y_pred.shape) < 3:
        cov = covariance_rowwise(x, y_pred)
        print(cov.shape)
    else:
        cov = np.zeros((x.shape[0], x.shape[1], x.shape[-1]))
        start_time = time.time()
        for i in range(x.shape[0]):
            print(i, time.time() - start_time)
            for j in range(x.shape[1]):
                cov[i, j, :] = covariance_rowwise(x[i, j, :, :],
                                                  y_pred[i, j, :])
    return cov


def build_df_for_plot(res, type_name):
    '''helper function to build dataframe for seaborn plot

    Args:
        res (ndarray): array for dataframe build
        type_name (list): name of each methods for plot

    Returns:
        dataframe: dataframe for seaborn plot
    '''
    base_dict = {}
    for j in range(len(type_name)):
        base_dict[type_name[j]] = res[j, :]
    df = pd.DataFrame(base_dict)
    mdf = pd.melt(df, var_name=['Type'])
    return mdf


def add_test_indicator(ax, lx, ly, rx, ry, data, dh, dx, ax_y0, ax_y1):
    '''plot predicted value

    Args:
        ax (matplotlib.pyplot.axis): axis for plot
        lx (float): x position of left side
        ly (float): max correlation (y axis) of left side
        rx (float): x position of right side
        ry (float): max correlation (y axis) of right side
        data (float): p value
        dh (float): offset to y axis
        dx (float): offset to x axis
        ax_y0 (float): lower limit of y axis
        ax_y1 (float): upper limit of y axis

    Returns:
        None
    '''
    flag_star = 1
    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05
        while data < p:
            text += '*'
            p /= 10.
        if p < 0.05:
            text = '*'
        if len(text) == 0:
            text = 'n.s.'
            flag_star = 0

    lx = dx + lx
    rx = dx + rx
    dh = dh * (ax_y1 - ax_y0)
    barh = 0.025 * (ax_y1 - ax_y0)
    y = max(ly, ry) + dh
    barx = [lx, rx]
    bary = [y, y]
    ax.plot(barx, bary, c='black')

    if flag_star:
        mid = ((lx + rx) / 2, y - barh * 1.2)
    else:
        mid = ((lx + rx) / 2, y - barh)

    kwargs = dict(ha='center', va='bottom', backgroundcolor='white')
    kwargs['fontsize'] = 24
    ax.text(*mid, text, **kwargs)


def plot_df(df, name, corr, p_value, dh, args, fontsize=24):
    '''plot boxplot for the PNF

    Args:
        df (dataframe): dataframe for seaborn plot
        name (str): name for saving
        corr (ndarray): correlation array for calculation
        p_value (ndarray): p value for several PNF
        dh (ndarray): offset for plot
        args: args from command line
        fontsize (optional, int): fontsize of plot

    Returns:
        None
    '''
    _, ax = plt.subplots(figsize=(15, 10), dpi=200)
    meanpointprops = dict(
        marker='^',
        markeredgecolor='black',
        markerfacecolor='white',
        markersize=20)
    width = 0.4
    colorp = ['#E3F2FD', '#64B5F6', '#1976D2', '#E91E63']
    colorp = sns.color_palette(colorp)
    sns.boxplot(
        ax=ax,
        x="Type",
        y="value",
        data=df,
        width=width,
        meanprops=meanpointprops,
        showmeans=True,
        palette=colorp)
    sns.swarmplot(ax=ax, x="Type", y="value", data=df, color=".25")
    sns.despine()
    dx = -0.2 + width / 2
    ax_y0, ax_y1 = ax.get_ylim()
    for i in [1, 2, 3]:
        if i == 1:
            rx = i - 0.05
        else:
            rx = i
        add_test_indicator(ax, 0, max(corr[0, :]), rx, max(corr[i, :]),
                           p_value[i - 1], dh[i - 1], dx, ax_y0, ax_y1)
    add_test_indicator(ax, 1, max(corr[0, :]), 2 - 0.05, max(corr[1, :]),
                       p_value[-1], dh[0], dx + 0.05, ax_y0, ax_y1)
    ax.set(xlabel=None)
    ax.set_ylabel('Agreement (Correlation) with Ground Truth')
    for item in ([ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() +
                 ax.get_yticklabels()):
        item.set_fontsize(fontsize)
        item.set_fontweight('medium')
    plt_name = os.path.join(args.out_dir, name + '_haufe_corr.png')
    plt.savefig(plt_name, bbox_inches='tight')


def run_test(corr):
    '''Perform t test for correlation

    Args:
        corr (ndarray): correlation array for calculation

    Returns:
        ndarray: p value calculated
    '''
    p = np.zeros((corr.shape[0], corr.shape[0]))
    for i in range(corr.shape[0]):
        for j in range(corr.shape[0]):
            if i == j:
                continue
            p[i, j] = ttest_rel(corr[i, :], corr[j, :])[1]
    print(p)
    return p


def plot_corr_box(corr, args):
    '''Plot correlation with boxplot

    Args:
        corr (ndarray): correlation array for plot
        args: args from command line

    Returns:
        None
    '''
    corr = np.swapaxes(corr[:, [3, 2, 0, 1], :], 0, 1)
    type_name = [
        'Basic\nMeta-matching\n(DNN) Training', 'Basic\nMeta-matching\n(DNN)',
        'Advanced\nMeta-matching\n(Stacking)', 'Classical (KRR)'
    ]
    corr = np.squeeze(np.abs(np.mean(corr, 2)))
    df = build_df_for_plot(corr, type_name)
    p_val = run_test(corr)
    plot_df(df, '35_phenotype' + args.stack_stem, corr,
            np.append(p_val[0, 1:], p_val[1, 2]), [0.085, 0.16, 0.235], args)


def haufe_transform_check(args):
    '''main function for Haufe transform check

    Args:
        args: args from command line

    Returns:
        None
    '''
    res_plot_npz = os.path.join(args.large_data_dir, 'haufe_res_for_plot.npz')
    if args.plot_only and os.path.isfile(res_plot_npz):
        npz = np.load(res_plot_npz)
        res_corr = npz['res_corr']
    else:
        print('\nCBIG meta-matching (DNN and DNN finetune) with argument: ' +
              str(args))

        # set all the seed
        seed = args.seed
        random.seed(seed)
        np.random.seed(seed)

        # load original data
        npz = os.path.join(args.in_dir, 'ukbb_dnn_input_cross_dataset.npz')
        npz = np.load(npz, allow_pickle=True)
        x_test = npz['x_test']
        tes_phe = npz['tes_phe']
        x_train = npz['x_train_raw']

        # npz save PNFs, for recover in case of fail due to long compute time
        res_npz_file = os.path.join(args.large_data_dir, 'haufe_res_PNF.npz')
        res_npz = {}
        if os.path.isfile(res_npz_file):
            res_npz = np.load(
                res_npz_file, allow_pickle=True)['arr_0'].tolist()
            match_phe = res_npz['match_phe'].astype(int)
        else:
            # load matched phenotype from DNN MM result
            npz = os.path.join(args.out_dir,
                               'meta_result_test_across_dataset.npz')
            npz = np.load(npz, allow_pickle=True)
            meta_phe = npz['meta_phe']
            meta_phe = np.squeeze(meta_phe[:, 3, :]).astype(int)  # k=100 only
            match_phe = np.zeros((meta_phe.shape[-1]))
            for i in range(meta_phe.shape[-1]):
                match_phe[i] = np.bincount(meta_phe[:, i]).argmax()
            match_phe = match_phe.astype(int)
            res_npz['match_phe'] = match_phe
            np.savez(res_npz_file, res_npz)

        # compute each PNF
        # PNF A: Basic Meta-matching (DNN) Training
        if 'pnf_A' in res_npz:
            pnf_A = res_npz['pnf_A']
        else:
            print("Compute PNF A")
            npz = os.path.join(args.out_dir, 'haufe_y_pred_train.npz')
            npz = np.load(npz)
            y_pred_train = npz['y_pred_train']
            x_train = npz['x_train']

            print(y_pred_train.shape, x_train.shape)
            pnf_A = compute_PNF(x_train, y_pred_train)
            res_npz['pnf_A'] = pnf_A
            np.savez(res_npz_file, res_npz)

        # PNF B: Basic Meta-matching (DNN)
        npz_stack_name = 'haufe_y_pred_100_stacking' + args.stack_stem + '.npz'
        if 'pnf_B' in res_npz:
            pnf_B = res_npz['pnf_B']
        else:
            print("Compute PNF B")
            npz = os.path.join(args.large_data_dir,
                               'haufe_y_pred_100_basic_mm.npz')
            npz = np.load(npz)
            basic_mm_pred_k_100 = npz['meta_pred_k_100']
            if 'stack_x_k_100' not in locals():
                npz = os.path.join(args.large_data_dir, npz_stack_name)
                npz = np.load(npz)
                stack_x_k_100 = npz['meta_x_k_100']
            pnf_B = compute_PNF(stack_x_k_100, basic_mm_pred_k_100)
            res_npz['pnf_B'] = pnf_B
            np.savez(res_npz_file, res_npz)

        # PNF C: Advanced Meta-matching (stacking)
        if 'pnf_C' in res_npz:
            pnf_C = res_npz['pnf_C']
        else:
            print("Compute PNF C")
            npz = os.path.join(args.large_data_dir, npz_stack_name)
            npz = np.load(npz)
            stack_pred_k_100 = npz['meta_pred_k_100']
            stack_x_k_100 = npz['meta_x_k_100']

            pnf_C = compute_PNF(stack_x_k_100, stack_pred_k_100)
            res_npz['pnf_C'] = pnf_C
            np.savez(res_npz_file, res_npz)

        # PNF D: "ground truth" KRR on full 1019 HCP dataset
        if 'pnf_D' in res_npz:
            pnf_D = res_npz['pnf_D']
        else:
            print("Compute PNF D")
            npz = os.path.join(args.inter_dir, 'HCP_haufe_all.npz')
            npz = np.load(npz, allow_pickle=True)
            pred_all_dict = npz['pred_all_dict'].item()
            HCP_krr_pred_1019 = []
            for i in tes_phe:
                HCP_krr_pred_1019.append(pred_all_dict[i])
            HCP_krr_pred_1019 = np.squeeze(HCP_krr_pred_1019)
            print(HCP_krr_pred_1019.shape)
            pnf_D = compute_PNF(x_test, np.swapaxes(HCP_krr_pred_1019, 0, 1))
            res_npz['pnf_D'] = pnf_D
            np.savez(res_npz_file, res_npz)

        # PNF E: classical KRR on 100 shot
        if 'pnf_E' in res_npz:
            pnf_E = res_npz['pnf_E']
        else:
            print("Compute PNF E")
            npz = os.path.join(args.inter_dir, 'HCP_haufe_100.npz')
            npz = np.load(npz, allow_pickle=True)
            pred_haufe_dict = npz['pred_haufe_dict'].item()
            HCP_krr_pred_100 = []
            for i in tes_phe:
                # print(pred_haufe_dict[i].shape)
                HCP_krr_pred_100.append(pred_haufe_dict[i])
            HCP_krr_pred_100 = np.squeeze(HCP_krr_pred_100)
            print(HCP_krr_pred_100.shape)
            if 'stack_x_k_100' not in locals():
                npz = os.path.join(args.large_data_dir, npz_stack_name)
                npz = np.load(npz)
                stack_x_k_100 = npz['meta_x_k_100']
            pnf_E = compute_PNF(stack_x_k_100,
                                np.swapaxes(HCP_krr_pred_100, 0, 1))
            res_npz['pnf_E'] = pnf_E
            np.savez(res_npz_file, res_npz)

        tmp_A = np.zeros((pnf_A.shape[0], len(tes_phe)))
        for i in range(len(tes_phe)):
            tmp_A[:, i] = pnf_A[:, match_phe[i]]

        print(match_phe.shape, pnf_A.shape, tmp_A.shape, pnf_B.shape,
              pnf_C.shape, pnf_D.shape, pnf_E.shape)

        res_corr = np.zeros((len(tes_phe), 4, 100))
        for i in range(len(tes_phe)):
            for j in range(100):
                res_corr[i, 0, j] += pearsonr(
                    np.squeeze(pnf_C[j, i, :]), pnf_D[:, i])[0]
                res_corr[i, 1, j] += pearsonr(
                    np.squeeze(pnf_E[j, i, :]), pnf_D[:, i])[0]
                res_corr[i, 2, j] += pearsonr(
                    np.squeeze(pnf_B[j, i, :]), pnf_D[:, i])[0]
                res_corr[i, 3, j] += pearsonr(tmp_A[:, i], pnf_D[:, i])[0]
            print(i, np.mean(res_corr[i, :, :], -1))

        print(np.mean(np.abs(np.mean(res_corr, -1)), 0))
        np.savez(res_plot_npz, res_corr=res_corr)

    plot_corr_box(res_corr, args)
    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--out_dir', '-o', type=str, default=config.OUT_DIR)
    parser.add_argument('--in_dir', type=str, default=config.IN_DIR)
    parser.add_argument('--inter_dir', type=str, default=config.INTER_DIR)
    parser.add_argument(
        '--large_data_dir', type=str, default=config.LARGE_DATA_DIR)
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--stack_stem', type=str, default='_restricted_alpha')
    parser.add_argument('--plot_only', type=bool, default=False)

    return parser.parse_args()


if __name__ == '__main__':
    haufe_transform_check(get_args())
