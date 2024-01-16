#!/usr/bin/env python
# coding: utf-8
'''
# Generate figures for the CBIG_MMP project
This python script contains scripts to generate figures in the paper
"Comparison of individualized behavioral predictions across anatomical,
diffusion and functional connectivity MRI, Ooi et al., 2022".

### Prerequisites
This code assumes that you have saved all regression results from kernel ridge
regression, linear ridge regression and elastic net using the scripts
`CBIG_MMP_<dataset>_collate_result_wrapper.m`. This would have resulted in a
.mat file saving the results of each single-feature-type model result for each
component score and the "grand average" (Average over all original scores,
58 in the HCP and 36 in the ABCD).

### Running the code
Create an output directory to store the images for the figures relating to HCP,
and a directory for images for ABCD. Modify the results directory and the
output directory for the images in part 3 of the code to suit your setup.
The code then generates all figures in the output directory.

Written by Leon Ooi and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

#######################################################
# Part 1
# Import libraries, set up plotting functions and set up directories
#######################################################

import scipy.io as scio
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt
from matplotlib import pyplot
from matplotlib.patches import PathPatch

#######################################################
# Part 2
# Set up results directory and output directory for figures
#######################################################
# IMPORTANT: Please change here to where the results are saved
MMP_dir = '/home/leon_ooi/storage/Multimodal_prediction_project/replication/'
HCP_results_dir = MMP_dir + 'HCP/collated_results/'
ABCD_results_dir = MMP_dir + 'ABCD/collated_results/'
HCP_output_dir = MMP_dir + 'HCP/collated_results/images/'
ABCD_output_dir = MMP_dir + 'ABCD/collated_results/images/'

#######################################################
# Part 3
# Define useful utilities
#######################################################


def adjust_box_widths(g, fac):
    """
    Adjust the widths of a seaborn-generated boxplot.

    Inputs:
    g: Figure
    fac: How much space to allocate between boxes
    """

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5 * (xmin + xmax)
                xhalf = 0.5 * (xmax - xmin)

                # setting new width of box
                xmin_new = xmid - fac * xhalf
                xmax_new = xmid + fac * xhalf
                verts_sub[np.equal(verts_sub[:, 0], xmin), 0] = xmin_new
                verts_sub[np.equal(verts_sub[:, 0], xmax), 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(np.array_equal(l.get_xdata(), [xmin, xmax])):
                        l.set_xdata([xmin_new, xmax_new])


def create_boxplot(size, font_size, acc_df, palette_sel, ylimits, metric,
                   outdir):
    """
    Main function for creating boxplot.

    Inputs:
    size: Tuple containing size of figure - width followed by height
    font_size: Font size for all labels
    acc_df: Dataframe of accuracy results
    palette_sel: Color palette to use
    ylimits: Vertical limits
    metric: Accuracy label on the y axis
    outdir: Directory to save image to
    """

    # create boxplot
    fig, ax = pyplot.subplots(figsize=size)
    sns.set()
    sns.set_style("ticks", {
        "xtick.major.size": font_size,
        "ytick.major.size": font_size
    })
    # plot
    bplot = sns.boxplot(
        x='Component',
        y=metric,
        hue='Feature',
        data=acc_df,
        fliersize=1,
        palette=palette_sel,
        linewidth=0.7)
    label_size = font_size + 1
    if font_size < 9:
        label_size = font_size + 3
    bplot.set_xticklabels(bplot.get_xmajorticklabels(), size=label_size - 1)
    adjust_box_widths(fig, 0.8)
    # modify axes
    ax.set(ylim=ylimits)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel(metric, fontsize=label_size)
    # modify legend
    handles, labels = ax.get_legend_handles_labels()
    ncol = 1
    if len(labels) > 6:
        ncol = 2
    ax.legend(
        handles=handles[0:],
        labels=labels[0:],
        fontsize=font_size,
        ncol=ncol,
        labelspacing=0.3)
    # save
    plt.pyplot.tight_layout()
    fig.savefig(outdir, dpi=450)
    plt.pyplot.close(fig)


def plot_HCP_boxplot(data, labels, palette_sel, font_size, ylimits, size,
                     metric, avg, outdir):
    """
    Plot boxplots for the HCP data.

    Inputs:
    data: Struct field from the results mat file contain prediction results
    labels: labels for the legend, same order as results
    palette_sel: Color palette to use
    font_size: Font size for all labels
    ylimits: Vertical limits
    size: Tuple containing size of figure - width followed by height
    metric: Accuracy label on the y axis
    avg: 1 to plot the grand average. 0 otherwise.
    outdir: Directory to save image to
    """

    # whether to plot average over all scores or not
    if avg == 1:
        scores = ['Cognition', 'Dissastisfaction', 'Emotion', 'Grand Average']
    else:
        scores = ['Cognition', 'Dissastisfaction', 'Emotion']
    n_scores = np.arange(len(scores))

    # append each score to dataframe and reshape for boxplot
    for n in range(0, len(labels)):
        if n == 0:
            acc_df = pd.DataFrame(data[0][0][n_scores, :, n].T, columns=scores)
            acc_df = acc_df.melt(
                value_vars=scores, var_name='Component', value_name=metric)
            acc_df['Feature'] = labels[n]
        else:
            df_data = pd.DataFrame(
                data[0][0][n_scores, :, n].T, columns=scores)
            df_data = df_data.melt(
                value_vars=scores, var_name='Component', value_name=metric)
            df_data['Feature'] = labels[n]
            acc_df = pd.concat([acc_df, df_data], ignore_index=True)
    # plot
    create_boxplot(size, font_size, acc_df, palette_sel, ylimits, metric,
                   outdir)


def plot_ABCD_boxplot(data, labels, palette_sel, font_size, ylimits, size,
                      metric, avg, outdir):
    """
    Plot boxplots for the ABCD data.

    Inputs:
    data: Struct field from the results mat file contain prediction results
    labels: labels for the legend, same order as results
    palette_sel: Color palette to use
    font_size: Font size for all labels
    ylimits: Vertical limits
    size: Tuple containing size of figure - width followed by height
    metric: Accuracy label on the y axis
    avg: 1 to plot the grand average. 0 otherwise.
    outdir: Directory to save image to
    """

    # whether to plot average over all scores or not
    if avg == 1:
        scores = ['Cognition', 'Personality', 'Mental Health', 'Grand Average']
    else:
        scores = ['Cognition', 'Personality', 'Mental Health']
    n_scores = np.arange(len(scores))

    # append each score to dataframe and reshape for boxplot
    for n in range(0, len(labels)):
        if n == 0:
            acc_df = pd.DataFrame(data[0][0][n_scores, :, n].T, columns=scores)
            acc_df = acc_df.melt(
                value_vars=scores, var_name='Component', value_name=metric)
            acc_df['Feature'] = labels[n]
        else:
            df_data = pd.DataFrame(
                data[0][0][n_scores, :, n].T, columns=scores)
            df_data = df_data.melt(
                value_vars=scores, var_name='Component', value_name=metric)
            df_data['Feature'] = labels[n]
            acc_df = pd.concat([acc_df, df_data], ignore_index=True)
    # plot
    create_boxplot(size, font_size, acc_df, palette_sel, ylimits, metric,
                   outdir)


def create_busy_boxplot(size, font_size, acc_df, palette_sel, ylimits, metric,
                        outdir):
    """
    Function for creating boxplot in which the x-labels need to be rotated.

    Inputs:
    size: Tuple containing size of figure - width followed by height
    font_size: Font size for all labels
    acc_df: Dataframe of accuracy results
    palette_sel: Color palette to use
    ylimits: Vertical limits
    metric: Accuracy label on the y axis
    outdir: Directory to save image to
    """

    # create boxplot
    fig, ax = pyplot.subplots(figsize=size)
    sns.set()
    sns.set_style("ticks", {
        "xtick.major.size": font_size,
        "ytick.major.size": font_size
    })
    # plot
    sns.boxplot(
        x='Component',
        y=metric,
        hue='Feature',
        data=acc_df,
        fliersize=1,
        palette=palette_sel,
        linewidth=0.7)
    # modify axes
    ax.set(ylim=ylimits)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    label_size = font_size + 1
    if font_size < 9:
        label_size = font_size + 3
    ax.set_xlabel('', fontsize=label_size)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha='right')
    ax.set_ylabel(metric, fontsize=label_size)
    # modify legend
    handles, labels = ax.get_legend_handles_labels()
    ncol = 1
    ax.legend(
        handles=handles[0:],
        labels=labels[0:],
        fontsize=font_size,
        ncol=ncol,
        labelspacing=0.3,
        bbox_to_anchor=(1.1, 1.1))
    # save
    plt.pyplot.tight_layout()
    fig.savefig(outdir, dpi=450)
    plt.pyplot.close(fig)


#######################################################
# Part 4
# This section generates plots from the main text.
#######################################################

# Figure 1: Mean KRR predictive performance from single modalities
# settings
font_size = 10
ylim = (-0.1, 0.6)
figsize = (6, 4)
avg = 1
metric = "Correlation"
labels = [
    'Anatomical (T1)', 'TBSS (dMRI)', 'Structural connectivity (dMRI)',
    'Functional connectivity (fMRI)'
]
palette = ["#95a5a6", "#9ACD32", "#04B8D2", "#e74c3c"]

# load HCP KRR results
krr_res = scio.loadmat(HCP_results_dir + 'KRR_corr_results.mat')
krr_mean = krr_res['results']['mean']
# mean
plot_HCP_boxplot(krr_mean, labels, palette, font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "KRR_matplotlib_mean.png")

# load ABCD KRR results
krr_res = scio.loadmat(ABCD_results_dir + 'KRR_corr_results.mat')
krr_mean = krr_res['results']['mean']
# mean
plot_ABCD_boxplot(krr_mean, labels, palette, font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "KRR_matplotlib_mean.png")

# Figure 2: Mean LRR predictive performance from single modalities
# settings
font_size = 10
ylim = (-0.1, 0.6)
figsize = (6, 4)
avg = 1
metric = "Correlation"
labels = [
    'Anatomical (T1)', 'TBSS (dMRI)', 'Structural connectivity (dMRI)',
    'Functional connectivity (fMRI)'
]
palette = ["#95a5a6", "#9ACD32", "#04B8D2", "#e74c3c"]

# load HCP LRR results
lrr_res = scio.loadmat(HCP_results_dir + 'LRR_corr_results.mat')
lrr_mean = lrr_res['results']['mean']
# mean
plot_HCP_boxplot(lrr_mean, labels, palette, font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "LRR_matplotlib_mean.png")

# load ABCD LRR results
lrr_res = scio.loadmat(ABCD_results_dir + 'LRR_corr_results.mat')
lrr_mean = lrr_res['results']['mean']
# mean
plot_ABCD_boxplot(lrr_mean, labels, palette, font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "LRR_matplotlib_mean.png")

# Figure 3: Mean elasticnet predictive performance from single modalities
# settings
font_size = 10
ylim = (-0.1, 0.6)
figsize = (6, 4)
avg = 1
metric = "Correlation"
labels = [
    'Anatomical (T1)', 'TBSS (dMRI)', 'Structural connectivity (dMRI)',
    'Functional connectivity (fMRI)'
]
palette = ["#95a5a6", "#9ACD32", "#04B8D2", "#e74c3c"]

# load HCP elasticnet results
en_res = scio.loadmat(HCP_results_dir + '/Elasticnet_corr_results.mat')
en_mean = en_res['results']['mean']
# mean
plot_HCP_boxplot(en_mean, labels, palette, font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "EN_matplotlib_mean.png")

# load ABCD elasticnet results
en_res = scio.loadmat(ABCD_results_dir + 'Elasticnet_corr_results.mat')
en_mean = en_res['results']['mean']
# mean
plot_ABCD_boxplot(en_mean, labels, palette, font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "EN_matplotlib_mean.png")

# Figure 4: Plot KRR best data
# settings
font_size = 10
ylim = (-0.1, 0.7)
figsize = (6, 4)
avg = 1
palette = ["#95a5a6", "#9ACD32", "#04B8D2", "#e74c3c"]
metric = "Correlation"
labels = [
    'Anatomical (T1)', 'TBSS (dMRI)', 'Structural connectivity (dMRI)',
    'Functional connectivity (fMRI)'
]

# load HCP KRR results
krr_res = scio.loadmat(HCP_results_dir + 'KRR_corr_results.mat')
krr_best = krr_res['results']['best']
# best
plot_HCP_boxplot(krr_best, labels, palette, font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "KRR_matplotlib_best.png")

# load ABCD KRR results
krr_res = scio.loadmat(ABCD_results_dir + 'KRR_corr_results.mat')
krr_best = krr_res['results']['best']
# best
plot_ABCD_boxplot(krr_best, labels, palette, font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "KRR_matplotlib_best.png")

# Figure 5: Plot HCP KRR data
# load HCP KRR results
krr_res = scio.loadmat(HCP_results_dir + 'KRR_corr_results.mat')
krr_t1 = krr_res['results']['t1']
krr_tbss = krr_res['results']['tbss']
krr_tractography = krr_res['results']['tractography']
krr_fmri = krr_res['results']['fmri']
krr_best = krr_res['results']['best']

# settings
font_size = 8
ylim = (-0.05, 0.7)
avg = 0
metric = "Correlation"
figsize = (4, 3.5)

# fMRI
labels = [
    'Resting', 'Social', 'Gambling', 'Language', 'Working Memory', 'Motor'
]
plot_HCP_boxplot(krr_fmri, labels, "OrRd", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "KRR_matplotlib_fmri.png")

# change y limits for other plots
ylim = (-0.1, 0.4)

# structural plot
labels = ['Thickness', 'Area', 'Volume']
plot_HCP_boxplot(krr_t1, labels, "Greys", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "KRR_matplotlib_struct.png")

# change fig size for remaining plots
figsize = (7.75, 3.5)
font_size = 10

# TBSS
labels = ['FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF']
plot_HCP_boxplot(krr_tbss, labels, "YlGn", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "KRR_matplotlib_tbss.png")

# MRtrix
labels = [
    'FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF', 'Stream count (log)',
    'Stream length'
]
plot_HCP_boxplot(krr_tractography, labels, "PuBu", font_size, ylim, figsize,
                 metric, avg, HCP_output_dir + "KRR_matplotlib_sc.png")

# Figure 6: Plot ABCD KRR data
# load ABCD KRR results
krr_res = scio.loadmat(ABCD_results_dir + 'KRR_corr_results.mat')
krr_t1 = krr_res['results']['t1']
krr_tbss = krr_res['results']['tbss']
krr_tractography = krr_res['results']['tractography']
krr_fmri = krr_res['results']['fmri']
krr_best = krr_res['results']['best']
krr_mean = krr_res['results']['mean']

# change fig size
# settings
font_size = 8
ylim = (-0.05, 0.7)
avg = 0
metric = "Correlation"
figsize = (4, 3.5)

# fMRI
labels = ['Resting', 'N-back', 'MID', 'SST']
plot_ABCD_boxplot(krr_fmri, labels, "OrRd", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "KRR_matplotlib_fmri.png")

# change y limits for other plots
ylim = (-0.15, 0.4)

# Structural
labels = ['Thickness', 'Area', 'Volume']
plot_ABCD_boxplot(krr_t1, labels, "Greys", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "KRR_matplotlib_struct.png")

# change fig size for remaining plots
figsize = (7.75, 3.5)
font_size = 10

# TBSS
labels = ['FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF']
plot_ABCD_boxplot(krr_tbss, labels, "YlGn", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "KRR_matplotlib_tbss.png")

# MRtrix
labels = [
    'FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF', 'Stream count (log)',
    'Stream length'
]
plot_ABCD_boxplot(krr_tractography, labels, "PuBu", font_size, ylim, figsize,
                  metric, avg, ABCD_output_dir + "KRR_matplotlib_sc.png")

# Figure 7: Plot HCP and ABCD combined models (Correlation)
# settings
font_size = 10
ylim = (-0.1, 0.75)
avg = 0
figsize = (6, 4)
palette = ["#FA8072", "#029386", "#ED0DD9", "#DBB40C"]
metric = "Correlation"

comb_res = scio.loadmat(HCP_results_dir + 'combined_models_corr_results.mat')
comb = comb_res['results']['combined']
# combined
labels = [
    'Best single-feature-type (KRR)', 'MultiKRR - All FC features',
    'Stacking - All FC models', 'Stacking - All MRI models'
]
plot_HCP_boxplot(comb, labels, palette, font_size, (-0.1, 0.8), figsize,
                 metric, avg,
                 HCP_output_dir + "KRR_matplotlib_combined_corr.png")

comb_res = scio.loadmat(ABCD_results_dir + 'combined_models_corr_results.mat')
comb = comb_res['results']['combined']
# combined
labels = [
    'Best single-feature-type (KRR)', 'MultiKRR - All FC features',
    'Stacking - All FC models', 'Stacking - All MRI models'
]
plot_ABCD_boxplot(comb, labels, palette, font_size, ylim, figsize, metric, avg,
                  ABCD_output_dir + "KRR_matplotlib_combined_corr.png")

#######################################################
# Part 5
# This section generates figures from the supplemental section.
#######################################################

# Figure S1: Mean KRR COD from single modalities
# settings
font_size = 10
ylim = (-0.1, 0.3)
figsize = (6, 4)
avg = 1
metric = "COD"
labels = [
    'Anatomical (T1)', 'TBSS (dMRI)', 'Structural connectivity (dMRI)',
    'Functional connectivity (fMRI)'
]
palette = ["#95a5a6", "#9ACD32", "#04B8D2", "#e74c3c"]

# load HCP KRR results
krr_res = scio.loadmat(HCP_results_dir + 'KRR_COD_results.mat')
krr_mean = krr_res['results']['mean']
# mean
plot_HCP_boxplot(krr_mean, labels, palette, font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "KRR_matplotlib_COD_mean.png")

# load ABCD KRR results
krr_res = scio.loadmat(ABCD_results_dir + 'KRR_COD_results.mat')
krr_mean = krr_res['results']['mean']
# mean
plot_ABCD_boxplot(krr_mean, labels, palette, font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "KRR_matplotlib_COD_mean.png")

# Figure S2-4: Mean HCP KRR (Every behaviour)
# Note: for these figures, please provide a list of real names of behaviors
# matching the order in the mat file.
# settings
font_size = 10
ylim = (-0.1, 0.7)
figsize = (7, 4)
metric = "Correlation"
labels = [
    'Anatomical (T1)', 'TBSS (dMRI)', 'Structural connectivity (dMRI)',
    'Functional connectivity (fMRI)'
]
palette = ["#95a5a6", "#9ACD32", "#04B8D2", "#e74c3c"]

# load HCP KRR results
krr_res = scio.loadmat(HCP_results_dir + 'KRR_corr_results.mat')
krr_mean_all = krr_res['results']['mean_all']

# use real names to plot: Please provide list of variable names
# matching the order in the mat file
scores = np.genfromtxt(
    'HCP_variables_to_predict_real_names.txt', dtype=str, delimiter='\n')

# append each score to dataframe and reshape for boxplot
for behav_set in range(0, 6):
    behav_max = (behav_set * 10) + 9
    behav_idx = range(behav_set * 10, min(behav_max, 58))
    print(behav_idx)
    subset_scores = scores[behav_idx]
    for n in range(0, len(labels)):
        if n == 0:
            acc_df = pd.DataFrame(
                krr_mean_all[0][0][behav_idx, :, n].T, columns=subset_scores)
            acc_df = acc_df.melt(
                value_vars=subset_scores,
                var_name='Component',
                value_name=metric)
            acc_df['Feature'] = labels[n]
        else:
            df_data = pd.DataFrame(
                krr_mean_all[0][0][behav_idx, :, n].T, columns=subset_scores)
            df_data = df_data.melt(
                value_vars=subset_scores,
                var_name='Component',
                value_name=metric)
            df_data['Feature'] = labels[n]
            acc_df = pd.concat([acc_df, df_data], ignore_index=True)

    create_busy_boxplot(
        figsize, font_size, acc_df, palette, ylim, metric,
        HCP_output_dir + "KRR_matplotlib_mean_all_" + str(behav_set) + ".png")

# Figure S5,6: Mean ABCD KRR (Every behaviour)
# Note: for these figures, please provide a list of real names of behaviors
# matching the order in the mat file.
# settings
font_size = 10
ylim = (-0.1, 0.7)
figsize = (7, 4)
metric = "Correlation"
labels = [
    'Anatomical (T1)', 'TBSS (dMRI)', 'Structural connectivity (dMRI)',
    'Functional connectivity (fMRI)'
]
palette = ["#95a5a6", "#9ACD32", "#04B8D2", "#e74c3c"]

# load HCP KRR results
krr_res = scio.loadmat(ABCD_results_dir + '/KRR_corr_results.mat')
krr_mean_all = krr_res['results']['mean_all']

# use real names to plot: Please provide list of variable names
# matching the order in the mat file
scores = np.genfromtxt(
    'ABCD_variables_to_predict_real_names.txt', dtype=str, delimiter='\n')

# append each score to dataframe and reshape for boxplot
for behav_set in range(0, 4):
    behav_max = (behav_set * 10) + 9
    behav_idx = range(behav_set * 10, min(behav_max, 36))
    print(behav_idx)
    subset_scores = scores[behav_idx]
    for n in range(0, len(labels)):
        if n == 0:
            acc_df = pd.DataFrame(
                krr_mean_all[0][0][behav_idx, :, n].T, columns=subset_scores)
            acc_df = acc_df.melt(
                value_vars=subset_scores,
                var_name='Component',
                value_name=metric)
            acc_df['Feature'] = labels[n]
        else:
            df_data = pd.DataFrame(
                krr_mean_all[0][0][behav_idx, :, n].T, columns=subset_scores)
            df_data = df_data.melt(
                value_vars=subset_scores,
                var_name='Component',
                value_name=metric)
            df_data['Feature'] = labels[n]
            acc_df = pd.concat([acc_df, df_data], ignore_index=True)

    create_busy_boxplot(
        figsize, font_size, acc_df, palette, ylim, metric,
        ABCD_output_dir + "KRR_matplotlib_mean_all_" + str(behav_set) + ".png")

# Figure S7: Plot LRR best data
# settings
font_size = 10
ylim = (-0.15, 0.65)
figsize = (6, 4)
avg = 1
palette = ["#95a5a6", "#9ACD32", "#04B8D2", "#e74c3c"]
metric = "Correlation"
labels = [
    'Anatomical (T1)', 'TBSS (dMRI)', 'Structural connectivity (dMRI)',
    'Functional connectivity (fMRI)'
]

# load HCP KRR results
lrr_res = scio.loadmat(HCP_results_dir + 'LRR_corr_results.mat')
lrr_best = lrr_res['results']['best']
# best
plot_HCP_boxplot(lrr_best, labels, palette, font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "LRR_matplotlib_best.png")

# load ABCD KRR results
lrr_res = scio.loadmat(ABCD_results_dir + 'LRR_corr_results.mat')
lrr_best = lrr_res['results']['best']
# best
plot_ABCD_boxplot(lrr_best, labels, palette, font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "LRR_matplotlib_best.png")

# Figure S8: Plot Elasticnet best data
# settings
font_size = 10
ylim = (-0.1, 0.7)
figsize = (6, 4)
avg = 1
palette = ["#95a5a6", "#9ACD32", "#04B8D2", "#e74c3c"]
metric = "Correlation"
labels = [
    'Anatomical (T1)', 'TBSS (dMRI)', 'Structural connectivity (dMRI)',
    'Functional connectivity (fMRI)'
]

# load HCP KRR results
en_res = scio.loadmat(HCP_results_dir + 'Elasticnet_corr_results.mat')
en_best = en_res['results']['best']
# best
plot_HCP_boxplot(en_best, labels, palette, font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "EN_matplotlib_best.png")

# load ABCD KRR results
en_res = scio.loadmat(ABCD_results_dir + 'Elasticnet_corr_results.mat')
en_best = en_res['results']['best']
# best
plot_ABCD_boxplot(en_best, labels, palette, font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "EN_matplotlib_best.png")

# Figure S9: HCP individual KRR models (COD)
# load HCP KRR results
krr_res = scio.loadmat(HCP_results_dir + 'KRR_COD_results.mat')
krr_t1 = krr_res['results']['t1']
krr_tbss = krr_res['results']['tbss']
krr_tractography = krr_res['results']['tractography']
krr_fmri = krr_res['results']['fmri']
krr_best = krr_res['results']['best']
krr_mean = krr_res['results']['mean']

# settings
font_size = 8
ylim = (-0.15, 0.4)
avg = 0
figsize = (4, 3.5)
metric = "COD"

# fMRI
labels = [
    'Resting', 'Social', 'Gambling', 'Language', 'Working Memory', 'Motor'
]
plot_HCP_boxplot(krr_fmri, labels, "OrRd", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "KRR_COD_matplotlib_fmri.png")

# change y limits for other plots
ylim = (-0.1, 0.2)

# structural plot
labels = ['Thickness', 'Area', 'Volume']
plot_HCP_boxplot(krr_t1, labels, "Greys", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "KRR_COD_matplotlib_struct.png")

# change fig size for remaining plots
figsize = (7.75, 3.5)
font_size = 10

# TBSS
labels = ['FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF']
plot_HCP_boxplot(krr_tbss, labels, "YlGn", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "KRR_COD_matplotlib_tbss.png")

# MRtrix
labels = [
    'FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF', 'Stream count (log)',
    'Stream length'
]
plot_HCP_boxplot(krr_tractography, labels, "PuBu", font_size, ylim, figsize,
                 metric, avg, HCP_output_dir + "KRR_COD_matplotlib_sc.png")

# Figure S10: ABCD individual KRR models (COD)
# load ABCD KRR results
krr_res = scio.loadmat(ABCD_results_dir + 'KRR_COD_results.mat')
krr_t1 = krr_res['results']['t1']
krr_tbss = krr_res['results']['tbss']
krr_tractography = krr_res['results']['tractography']
krr_fmri = krr_res['results']['fmri']
krr_best = krr_res['results']['best']
krr_mean = krr_res['results']['mean']

# change fig size
# settings
font_size = 8
ylim = (-0.15, 0.4)
avg = 0
metric = "COD"
figsize = (4, 3.5)

# fMRI
labels = ['Resting', 'N-back', 'MID', 'SST']
plot_ABCD_boxplot(krr_fmri, labels, "OrRd", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "KRR_COD_matplotlib_fmri.png")

# change y limits for other plots
ylim = (-0.1, 0.2)

# Structural
labels = ['Thickness', 'Area', 'Volume']
plot_ABCD_boxplot(krr_t1, labels, "Greys", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "KRR_COD_matplotlib_struct.png")

# change fig size for remaining plots
figsize = (7.75, 3.5)
font_size = 10

# TBSS
labels = ['FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF']
plot_ABCD_boxplot(krr_tbss, labels, "YlGn", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "KRR_COD_matplotlib_tbss.png")

# MRtrix
labels = [
    'FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF', 'Stream count (log)',
    'Stream length'
]
plot_ABCD_boxplot(krr_tractography, labels, "PuBu", font_size, ylim, figsize,
                  metric, avg, ABCD_output_dir + "KRR_COD_matplotlib_sc.png")

# Figure S11: HCP individual LRR models
# load HCP LRR results
lrr_res = scio.loadmat(HCP_results_dir + 'LRR_corr_results.mat')
lrr_t1 = lrr_res['results']['t1']
lrr_tbss = lrr_res['results']['tbss']
lrr_tractography = lrr_res['results']['tractography']
lrr_fmri = lrr_res['results']['fmri']
lrr_best = lrr_res['results']['best']
lrr_mean = lrr_res['results']['mean']

# change fig size
# settings
font_size = 8
ylim = (-0.05, 0.7)
avg = 0
metric = "Correlation"
figsize = (4, 3.5)

# fMRI
labels = [
    'Resting', 'Social', 'Gambling', 'Language', 'Working Memory', 'Motor'
]
plot_HCP_boxplot(lrr_fmri, labels, "OrRd", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "LRR_matplotlib_fmri.png")

# change y limits for other plots
ylim = (-0.2, 0.4)

# structural plot
labels = ['Thickness', 'Area', 'Volume']
plot_HCP_boxplot(lrr_t1, labels, "Greys", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "LRR_matplotlib_struct.png")

# change fig size for remaining plots
figsize = (7.75, 3.5)
font_size = 10

# TBSS
labels = ['FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF']
plot_HCP_boxplot(lrr_tbss, labels, "YlGn", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "LRR_matplotlib_tbss.png")

# MRtrix
labels = [
    'FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF', 'Stream count (log)',
    'Stream length'
]
plot_HCP_boxplot(lrr_tractography, labels, "PuBu", font_size, ylim, figsize,
                 metric, avg, HCP_output_dir + "LRR_matplotlib_sc.png")

# Figure S12: HCP individual Elasticnet models
# load HCP elasticnet results
en_res = scio.loadmat(HCP_results_dir + 'Elasticnet_corr_results.mat')
en_t1 = en_res['results']['t1']
en_tbss = en_res['results']['tbss']
en_tractography = en_res['results']['tractography']
en_fmri = en_res['results']['fmri']
en_best = en_res['results']['best']
en_mean = en_res['results']['mean']

# change fig size
# settings
font_size = 8
ylim = (-0.05, 0.7)
avg = 0
metric = "Correlation"
figsize = (4, 3.5)

# fMRI
labels = [
    'Resting', 'Social', 'Gambling', 'Language', 'Working Memory', 'Motor'
]
plot_HCP_boxplot(en_fmri, labels, "OrRd", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "EN_matplotlib_fmri.png")

# change y limits for other plots
ylim = (-0.25, 0.4)

# structural plot
labels = ['Thickness', 'Area', 'Volume']
plot_HCP_boxplot(en_t1, labels, "Greys", font_size, ylim, figsize, metric, avg,
                 HCP_output_dir + "EN_matplotlib_struct.png")

# change fig size for remaining plots
figsize = (7.75, 3.5)
font_size = 10

# TBSS
labels = ['FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF']
plot_HCP_boxplot(en_tbss, labels, "YlGn", font_size, ylim, figsize, metric,
                 avg, HCP_output_dir + "EN_matplotlib_tbss.png")

# MRtrix
labels = [
    'FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF', 'Stream count (log)',
    'Stream length'
]
plot_HCP_boxplot(en_tractography, labels, "PuBu", font_size, ylim, figsize,
                 metric, avg, HCP_output_dir + "EN_matplotlib_sc.png")

# Figure S13: ABCD individual LRR models
# load HCP LRR results
lrr_res = scio.loadmat(ABCD_results_dir + 'LRR_corr_results.mat')
lrr_t1 = lrr_res['results']['t1']
lrr_tbss = lrr_res['results']['tbss']
lrr_tractography = lrr_res['results']['tractography']
lrr_fmri = lrr_res['results']['fmri']
lrr_best = lrr_res['results']['best']
lrr_mean = lrr_res['results']['mean']

# change fig size
# settings
font_size = 8
ylim = (-0.1, 0.6)
avg = 0
metric = "Correlation"
figsize = (4, 3.5)

# fMRI
labels = ['Resting', 'N-back', 'MID', 'SST']
plot_ABCD_boxplot(lrr_fmri, labels, "OrRd", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "LRR_matplotlib_fmri.png")

# change y limits for other plots
ylim = (-0.2, 0.4)
font_size = 10

# Structural
labels = ['Thickness', 'Area', 'Volume']
plot_ABCD_boxplot(lrr_t1, labels, "Greys", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "LRR_matplotlib_struct.png")

# change fig size for remaining plots
figsize = (7.75, 3.5)

# TBSS
labels = ['FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF']
plot_ABCD_boxplot(lrr_tbss, labels, "YlGn", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "LRR_matplotlib_tbss.png")

# MRtrix
labels = [
    'FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF', 'Stream count (log)',
    'Stream length'
]
plot_ABCD_boxplot(lrr_tractography, labels, "PuBu", font_size, ylim, figsize,
                  metric, avg, ABCD_output_dir + "LRR_matplotlib_sc.png")

# Figure S14: ABCD individual Elasticnet models
# load HCP elasticnet results
en_res = scio.loadmat(ABCD_results_dir + 'Elasticnet_corr_results.mat')
en_t1 = en_res['results']['t1']
en_tbss = en_res['results']['tbss']
en_tractography = en_res['results']['tractography']
en_fmri = en_res['results']['fmri']
en_best = en_res['results']['best']
en_mean = en_res['results']['mean']

# change fig size
# settings
font_size = 8
ylim = (-0.1, 0.6)
avg = 0
metric = "Correlation"
figsize = (4, 3.5)

# fMRI
labels = ['Resting', 'N-back', 'MID', 'SST']
plot_ABCD_boxplot(en_fmri, labels, "OrRd", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "EN_matplotlib_fmri.png")

# change y limits for other plots
ylim = (-0.2, 0.4)
font_size = 10

# Structural
labels = ['Thickness', 'Area', 'Volume']
plot_ABCD_boxplot(en_t1, labels, "Greys", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "EN_matplotlib_struct.png")

# change fig size for remaining plots
figsize = (7.75, 3.5)

# TBSS
labels = ['FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF']
plot_ABCD_boxplot(en_tbss, labels, "YlGn", font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "EN_matplotlib_tbss.png")

# MRtrix
labels = [
    'FA', 'MD', 'AD', 'RD', 'OD', 'ISOVF', 'ICVF', 'Stream count (log)',
    'Stream length'
]
plot_ABCD_boxplot(en_tractography, labels, "PuBu", font_size, ylim, figsize,
                  metric, avg, ABCD_output_dir + "EN_matplotlib_sc.png")

# Figure S15: Stacking best of each modality
# settings
font_size = 10
ylim = (-0.1, 0.75)
avg = 0
figsize = (6, 4)
palette = ["#FA8072", "#029386", "#ED0DD9", "#DBB40C"]
metric = "Correlation"

comb_res = scio.loadmat(HCP_results_dir + 'best_models_corr_results.mat')
comb = comb_res['results']['best']
# pick best features
recombine = np.zeros((3, 60, 4))
recombine[:, :, 0] = comb[0][0][:, :, 0]
recombine[0, :, 1] = comb[0][0][0, :, 2]
recombine[1, :, 1] = comb[0][0][1, :, 3]
recombine[2, :, 1] = comb[0][0][2, :, 4]
recombine[:, :, 2] = comb[0][0][:, :, 1]
recombine[:, :, 3] = comb[0][0][:, :, 5]
recombine = [[recombine]]
# combined
labels = [
    'Best single-feature-type (KRR)', 'Stacking - Best of each modality',
    'Stacking - All FC models', 'Stacking - All MRI models'
]
plot_HCP_boxplot(recombine, labels, palette, font_size, (-0.1, 0.8), figsize,
                 metric, avg,
                 HCP_output_dir + "stacking_matplotlib_best_corr.png")

comb_res = scio.loadmat(ABCD_results_dir + 'best_models_corr_results.mat')
comb = comb_res['results']['best']
# pick best features
recombine = np.zeros((3, 120, 4))
recombine[:, :, 0] = comb[0][0][:, :, 0]
recombine[0, :, 1] = comb[0][0][0, :, 2]
recombine[1, :, 1] = comb[0][0][1, :, 3]
recombine[2, :, 1] = comb[0][0][2, :, 4]
recombine[:, :, 2] = comb[0][0][:, :, 1]
recombine[:, :, 3] = comb[0][0][:, :, 5]
recombine = [[recombine]]
# combined
labels = [
    'Best single-feature-type (KRR)', 'Stacking - Best of each modality',
    'Stacking - All FC models', 'Stacking - All MRI models'
]
plot_ABCD_boxplot(recombine, labels, palette, font_size, ylim, figsize, metric,
                  avg, ABCD_output_dir + "stacking_matplotlib_best_corr.png")

# Figure S16: Plot HCP and ABCD combined models (COD)
# settings
font_size = 10
ylim = (-0.15, 0.55)
avg = 0
figsize = (6, 4)
palette = ["#FA8072", "#029386", "#ED0DD9", "#DBB40C"]
metric = "COD"

comb_res = scio.loadmat(HCP_results_dir + 'combined_models_COD_results.mat')
comb = comb_res['results']['combined']
# combined
labels = [
    'Best single-feature-type (KRR)', 'MultiKRR - All FC features',
    'Stacking - All FC models', 'Stacking - All MRI models'
]
plot_HCP_boxplot(comb, labels, palette, font_size, ylim, figsize, metric, avg,
                 HCP_output_dir + "KRR_matplotlib_combined_COD.png")

comb_res = scio.loadmat(ABCD_results_dir + 'combined_models_COD_results.mat')
comb = comb_res['results']['combined']
# combined
labels = [
    'Best single-feature-type (KRR)', 'MultiKRR - All FC features',
    'Stacking - All FC models', 'Stacking - All MRI models'
]
plot_ABCD_boxplot(comb, labels, palette, font_size, ylim, figsize, metric, avg,
                  ABCD_output_dir + "KRR_matplotlib_combined_COD.png")
