# coding: utf-8

# this script will perform a hierarchical clustering on the predictive network features

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


def cluster_mat(array,
                list_var,
                pdist_metric,
                link_metric,
                c_dict,
                thresh=3,
                path_out=None,
                fname=None):
    """
    This function is to perform a clustering and plot a correlation matrix and dendrogram.
    Args:
        array:          Input data array. Each column is an variable and each row is an oberservation. 
                     
        list_var:       List of the variable names in array. Must have the same length as the number of
                        columns in array.
        
        pdist_metric:   Metric used to compute pairwise distance of variables. see sch.distance.pdist function
                        for acceptable values.
                        
        link_metric:    See sch.linkage function for acceptable values for method and metric.
        
        c_dict:         Color dictionary to color the tick lables.
        
        thresh:         Threshold used to cut off the clusters based on distance. Default is 3.
        
        path_out:       Output directory to save the figures. If None, No figure will be saved.
        
        fname:          File name for the saved figures.
    Returns:
        None
    """
    df = pd.DataFrame(data=array, columns=list_var)
    df_corr = df[list_var].corr()
    df_h = df_corr.copy()
    corr = df_h.values
    # perform the clustering
    # see sch.distance.pdist function for acceptable values for metric
    pdist = sch.distance.pdist(corr, metric=pdist_metric)
    # see sch.linkage function for acceptable values for method and metric
    linkage = sch.linkage(pdist, method=link_metric, metric=pdist_metric)
    idx = sch.fcluster(linkage, t=3, criterion='maxclust')

    # plot dendrogram
    plt.figure(figsize=(12, 4))
    plt.ylabel('distance')
    r = sch.dendrogram(
        linkage, color_threshold=thresh, labels=list_var, leaf_rotation=90)
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        lbl.set_color(c_dict[lbl.get_text()])
    if isinstance(path_out, str):
        plt.savefig(
            path_out + '/dendrogram_' + fname + '.png', bbox_inches='tight')
        plt.savefig(
            path_out + '/dendrogram_' + fname + '.svg', bbox_inches='tight')

    # plot rearranged clustered correlation matrix
    cols = r['ivl']
    df_h = df_h.reindex_axis(cols, axis=1)
    df_h = df_h.reindex_axis(cols, axis=0)
    #plt.figure(figsize=(25,18.75))
    plt.figure(figsize=(10, 8))
    g = sns.heatmap(
        df_h,
        vmin=-0.8,
        vmax=0.8,
        cmap="BrBG",
        xticklabels=df_h.columns,
        yticklabels=df_h.columns,
        linewidth=0.5,
        linecolor='black')
    for tick_label in g.get_yticklabels():
        tick_text = tick_label.get_text()
        tick_label.set_color(c_dict[tick_text])
    for tick_label in g.get_xticklabels():
        tick_text = tick_label.get_text()
        tick_label.set_color(c_dict[tick_text])
    if isinstance(path_out, str):
        plt.savefig(
            path_out + '/sim_matrix_' + fname + '.png', bbox_inches='tight')
        plt.savefig(
            path_out + '/sim_matrix_' + fname + '.svg', bbox_inches='tight')
    return idx, df_h, pdist, cols


# function to plot a correlation matrix without clustering
def plot_raw_matrix(array, list_var, c_dict, size, path_out, fname):
    """
    This function is to plot a correlation matrix among variables.
    Args:
        array:          Input data array. Each column is an variable and each row is an oberservation. 
                     
        list_var:       List of the variable names in array. Must have the same length as the number of
                        columns in array.
        
        c_dict:         Color dictionary to color the tick lables.
        
        size:           Figure size.
        
        path_out:       Output directory to save the figures. If None, No figure will be saved.
        
        fname:          File name for the saved figures.
    Returns:
        None
    """
    df = pd.DataFrame(data=array, columns=list_var)
    df_corr = df[list_var].corr()
    plt.figure(figsize=size)
    g = sns.heatmap(
        df_corr,
        vmin=-0.8,
        vmax=0.8,
        cmap="BrBG",
        xticklabels=df_corr.columns,
        yticklabels=df_corr.columns,
        linewidth=0.3,
        linecolor='black')
    for tick_label in g.get_yticklabels():
        tick_text = tick_label.get_text()
        tick_label.set_color(c_dict[tick_text])
    for tick_label in g.get_xticklabels():
        tick_text = tick_label.get_text()
        tick_label.set_color(c_dict[tick_text])
    if isinstance(path_out, str):
        plt.savefig(
            path_out + '/sim_matrix_' + fname + '.png', bbox_inches='tight')
        plt.savefig(
            path_out + '/sim_matrix_' + fname + '.svg', bbox_inches='tight')


def rearrange_matrix_sig(array, idx):
    nn = 0
    new_mat = np.zeros((array.shape[0], len(idx)))
    for i in idx:
        new_mat[:, nn] = array[:, i]
        nn += 1
    return new_mat


# dictionary for colours of different behavioural scales
c_dict = {
    'vocabulary': 'red',
    'attention': 'red',
    'working memory': 'red',
    'executive function': 'red',
    'processing speed': 'red',
    'episodic memory': 'red',
    'reading': 'red',
    'fluid cognition': 'red',
    'crystallized cognition': 'red',
    'overall cognition': 'red',
    'short delay recall': 'firebrick',
    'long delay recall': 'firebrick',
    'fluid intelligence': 'orangered',
    'visuospatial accuracy': 'lightcoral',
    'visuospatial reaction time': 'lightcoral',
    'visuospatial efficiency': 'lightcoral',
    'anxious depressed': 'blue',
    'withdrawn depressed': 'blue',
    'somatic complaints': 'blue',
    'social problems': 'blue',
    'thought problems': 'blue',
    'attention problems': 'blue',
    'rulebreaking behavior': 'blue',
    'aggressive behavior': 'blue',
    'total psychosis symptoms': 'cornflowerblue',
    'psychosis severity': 'cornflowerblue',
    'mania severity': 'dodgerblue',
    'negative urgency': 'black',
    'lack of planning': 'black',
    'sensation seeking': 'black',
    'positive urgency': 'black',
    'lack perseverance': 'black',
    'behavioral inhibition': 'dimgrey',
    'reward responsiveness': 'dimgrey',
    'drive': 'dimgrey',
    'fun seeking': 'dimgrey'
}


def plot_PFM_clustering(in_dir, data_dir, score_predicted, out_dir):
    """
    This function is a wrapper function to run the clustering of predictive feature matrix (PFM).
    Args:
        in_dir:             Input directory that contains the list of all variables.
                     
        data_dir:           Data directory that contains the PFM files. This directory should
                            be the output directory of ../matrix_plots/CBIG_TRBPC_plot_avg_relevance.m
        
        score_predicted:    List of varaibles that are significantly predicted. Clustering will only be 
                            performed on the significantly predicted behaviors.
                        
        out_dir:            Output directory to store the output figures.

    Returns:
        None
    """
    # list of all behaviors
    with open(in_dir + '/variables_to_predict.txt') as file:
        list_var = file.read().splitlines()

    # list of behaviors that were significantly predicted
    with open(score_predicted) as file:
        list_sig_r = file.read().splitlines()

    idx_sig = [list_var.index(i) for i in list_sig_r]
    # plot data-driven clustering when all the fmri conditions are stacked together
    mat_stack = sio.loadmat(data_dir + '/stacked_relevance_vectors.mat')
    stacked_mat = rearrange_matrix_sig(mat_stack['stack'], idx_sig)
    [stacked_idx, stacked_df_h, pdist, sig_var_reorder] = cluster_mat(
        stacked_mat, list_sig_r, 'euclidean', 'average', c_dict, 2.5, out_dir,
        'datadriven_stacked')
    # plot hypothesis driven clustering
    plot_raw_matrix(stacked_mat, list_sig_r, c_dict, (10, 8), out_dir,
                    'hypothesis_stacked')
    return list_var, list_sig_r, sig_var_reorder


def plot_PFM_clustering_separate_condition(list_var, var_reorder, data_dir,
                                           out_dir, fname):
    """ plot the giant matrix of PFM similarity across all behaviors and conditions """
    # get the indices of cols from list_var
    idx_cols = [list_var.index(i) for i in var_reorder]

    # load data
    mat_fc = sio.loadmat(data_dir + '/relevance_vectors.mat')

    # rearrange each fMRI matrix according to the clustering from the stacked matrix (all 4 fMRI vectors together)
    rs_mat_s = rearrange_matrix_sig(mat_fc['struct_fc_vec']['rs'][0][0],
                                    idx_cols)
    mid_mat_s = rearrange_matrix_sig(mat_fc['struct_fc_vec']['mid'][0][0],
                                     idx_cols)
    nback_mat_s = rearrange_matrix_sig(mat_fc['struct_fc_vec']['nback'][0][0],
                                       idx_cols)
    sst_mat_s = rearrange_matrix_sig(mat_fc['struct_fc_vec']['sst'][0][0],
                                     idx_cols)

    # create np array with all the vectors ordered by behavior and then fmri condition
    concat_mat_clus = np.zeros((rs_mat_s.shape[0], rs_mat_s.shape[1] * 4))
    nn = 0
    for bb in range(rs_mat_s.shape[1]):  # for each sig. behavior
        concat_mat_clus[:, nn] = rs_mat_s[:, bb]  # put in rs
        nn += 1
        concat_mat_clus[:, nn] = mid_mat_s[:, bb]  # put in mid
        nn += 1
        concat_mat_clus[:, nn] = sst_mat_s[:, bb]  # put in sst
        nn += 1
        concat_mat_clus[:, nn] = nback_mat_s[:, bb]  # put in nback
        nn += 1

    # set up dictionary for coloring labels
    list_fmri = ['REST ', 'MID ', 'SST ', 'NBACK ']
    list_sig_r_fmri_clus = []
    for behav in var_reorder:
        for fmri in list_fmri:
            tmp = fmri + behav
            list_sig_r_fmri_clus.append(tmp)
    colours = []
    for behav in var_reorder:
        colours.append(c_dict[behav])
    colours_fmri = []
    for cc in colours:
        colours_fmri += [cc] * len(list_fmri)
    fmri_c_dict = dict(zip(list_sig_r_fmri_clus, colours_fmri))

    # plot and save
    plt.rcParams.update({'font.size': 10})
    plot_raw_matrix(concat_mat_clus, list_sig_r_fmri_clus, fmri_c_dict,
                    (25, 18.75), out_dir, fname)


if __name__ == "__main__":
    np.random.seed(seed=1)
    if len(sys.argv) < 5:
        raise Exception('not enough input arguments')

    in_dir = sys.argv[1]
    out_dir = sys.argv[2]
    data_dir = sys.argv[3]
    score_predicted = sys.argv[4]
    list_var, sig_var, sig_var_reorder = plot_PFM_clustering(
        in_dir, data_dir, score_predicted, out_dir)
    plot_PFM_clustering_separate_condition(list_var, sig_var_reorder, data_dir,
                                           out_dir, 'datadriven_ind')
    plot_PFM_clustering_separate_condition(list_var, sig_var, data_dir,
                                           out_dir, 'hypothesis_ind')
