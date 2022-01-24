# coding: utf-8

# # this notebook produces matrices showing the correlations among fmri conditions and behavioral clusters

# Written by Angela Tam & CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy import stats
import scipy.cluster.hierarchy as sch
import os
import sys


def plot_PFM_similarity_of_conditions(in_dir, out_dir, cluster):
    """
    This function is to plot the similarity of predictive feature matrix (PFM)
    of different brain states.
    Args:
        in_dir:      Input directory that contains the PFM files. This directory should
                     be the output directory of ../matrix_plots/CBIG_TRBPC_plot_avg_relevance.m
                     
        out_dir:     Output directory to store the output figures.
        
        cluster:     Type of the cluster. Choose from 'hypothesis' or 'datadriven'
    Returns:
        None
    """
    # load data
    # this .mat file is an output from ../matrix_plots/CBIG_TRBPC_plot_avg_relevance.m
    if cluster == "hypothesis":
        mat_fc = sio.loadmat(in_dir + '/datadriven/mean_vec_clus_fmri.mat')
    elif cluster == "datadriven":
        mat_fc = sio.loadmat(in_dir + '/hypothesis/mean_vec_clus_fmri.mat')
    else:
        print("cluster must be 'hypothesis' or 'datadriven'")
        return

    # rearrange the matrix so that the order is rest, MID, SST, Nback (swap SST and Nback)
    mat_r = np.zeros((87571, 12))
    mat_r = np.copy(mat_fc['vec_stacked'])
    mat_r[:, 2] = mat_fc['vec_stacked'][:, 3]
    mat_r[:, 3] = mat_fc['vec_stacked'][:, 2]
    mat_r[:, 6] = mat_fc['vec_stacked'][:, 7]
    mat_r[:, 7] = mat_fc['vec_stacked'][:, 6]
    mat_r[:, 10] = mat_fc['vec_stacked'][:, 11]
    mat_r[:, 11] = mat_fc['vec_stacked'][:, 10]

    # list of conditions
    list_col = [
        'Cognitive REST', 'Cognitive MID', 'Cognitive SST', 'Cognitive NBACK',
        'Personality REST', 'Personality MID', 'Personality SST',
        'Personality NBACK', 'Mental health REST', 'Mental health MID',
        'Mental health SST', 'Mental health NBACK'
    ]

    # build data frame and calculate correlations
    df_corr_clus = pd.DataFrame(data=mat_r, columns=list_col)
    df_corr_clus = df_corr_clus.corr()

    # plot correlations
    plt.figure(figsize=(10, 8))
    g = sns.heatmap(
        df_corr_clus,
        vmin=-0.8,
        vmax=0.8,
        cmap="BrBG",
        xticklabels=df_corr_clus.columns,
        yticklabels=df_corr_clus.columns,
        linewidth=0.5,
        linecolor='black')
    plt.savefig(out_dir + '/' + cluster + '_behav_fmri_sim_matrix.png')
    plt.savefig(out_dir + '/' + cluster + '_behav_fmri_sim_matrix.svg')


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise Exception('not enough input arguments')

    in_dir = sys.argv[1]
    out_dir = sys.argv[2]
    plot_PFM_similarity_of_conditions(in_dir, out_dir, 'datadriven')
    plot_PFM_similarity_of_conditions(in_dir, out_dir, 'hypothesis')
