# coding: utf-8

# # this function will perform a hierarchical clustering on the raw behavioral scores

# Written by Angela Tam & CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import scipy.cluster.hierarchy as sch
import os
import sys


def cluster_mat(df, list_var, c_dict, path_out, fname, thresh=3):
    """
    This function is to perform a clustering and plot a correlation matrix and dendrogram.
    Args:
        df:             Input data frame.
                     
        list_var:       List of the variable names in array. Must have the same length as the number of
                        columns in array.
        
        c_dict:         Color dictionary to color the tick lables.
        
        path_out:       Output directory to save the figures. If None, No figure will be saved.
        
        fname:          File name for the saved figures.
        
        thresh:         Threshold used to cut off the clusters based on distance. Default is 3.
    Returns:
        None
    """
    df_corr = df[list_var].corr()
    df_h = df_corr.copy()
    corr = df_h.values
    # perform the clustering
    pdist = sch.distance.pdist(corr)
    linkage = sch.linkage(pdist, method='average')
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
        plt.savefig(path_out + '/dendrogram_' + fname + '.svg')

    # plot rearranged clustered correlation matrix
    cols = r['ivl']
    df_h = df_h.reindex(cols, axis=1)
    df_h = df_h.reindex(cols, axis=0)
    plt.figure(figsize=(10, 8))
    g = sns.heatmap(
        df_h,
        vmin=-0.8,
        vmax=0.8,
        cmap="BrBG",
        xticklabels=df_h.columns,
        yticklabels=df_h.columns)
    for tick_label in g.get_yticklabels():
        tick_text = tick_label.get_text()
        tick_label.set_color(c_dict[tick_text])
    for tick_label in g.get_xticklabels():
        tick_text = tick_label.get_text()
        tick_label.set_color(c_dict[tick_text])
    if isinstance(path_out, str):
        plt.savefig(path_out + '/sim_matrix_' + fname + '.svg')
    return idx, cols, df_h


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
    'social problems': 'blue',
    'thought problems': 'blue',
    'attention problems': 'blue',
    'total psychosis symptoms': 'cornflowerblue',
    'psychosis severity': 'cornflowerblue',
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


def plot_similarity_of_behavior(sig_score, path_out):
    """
    This servers as a wrapper function to plot the similarity of ABCD raw behaviors.
    Args:
        sig_score:      Behaviors that are significantly predicted.

        path_out:       Output directory to save the figures.
    Returns:
        None
    """
    # set random seed
    np.random.seed(seed=1)

    # input directories
    # private input
    in_dir_priv = os.getenv(
        'CBIG_REPDATA_DIR'
    ) + '/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/figures/input'
    # public input
    in_dir_pub = os.getenv(
        'CBIG_CODE_DIR') + '/stable_projects/predict_phenotypes'
    in_dir_pub += '/ChenTam2022_TRBPC/figure_utilities/input'

    df_demog = pd.read_csv(
        in_dir_priv +
        '/abcd_2.0_demog_cog_neuropsych_completecases_groupsites.csv')

    # load list of behaviors (named by ABCD)
    with open(in_dir_pub + '/variables_to_predict_abcd.txt') as file:
        list_var = file.read().splitlines()

    # load list of behaviors (named by us in 'lay people' terms)
    with open(in_dir_pub + '/variables_to_predict.txt') as file:
        list_var_lay = file.read().splitlines()

    # load significantly predicted behaviors
    with open(sig_score) as file:
        list_sig = file.read().splitlines()

    # change the test names to more laypeople terms
    df_demog.rename(columns=dict(zip(list_var, list_var_lay)), inplace=True)

    # perform a hierarchical clustering
    [idx, list_r, df_corr_hi] = cluster_mat(
        df_demog[list_sig], list_sig, c_dict, path_out, 'behavior', thresh=1.6)


if __name__ == "__main__":

    if len(sys.argv) < 3:
        raise Exception('not enough input arguments')

    sig_score = sys.argv[1]
    out_dir = sys.argv[2]
    plot_similarity_of_behavior(sig_score, out_dir)
