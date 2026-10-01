'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
#############################################################################
# ===========================================================================
# IMPORTS
# ===========================================================================
#############################################################################
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import copy
import pickle

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
sys.path.append(os.path.join(PROJECT_DIR, 'util'))
import CBIG_TMSmotion_util as util


#############################################################################
# ===========================================================================
# WRAPPER FUNCTIONS
# ===========================================================================
#############################################################################
# Wrapper 1
def select_subID_sim(main_config_file, pl_config_file, dataset, subID):
    """Wrapper that filters whether a subID satisfies simulation-based selection requirements.

    Overall, this wrapper function filters which subIDs successfully
    fulfills necessary simulation requirements. For example, being able
    to generate the minimum # of simulation runs for each motion bin
    (num_sim_per_bin), and having >= minimum # of "sufficient" bins
    (num_sample_bin2select_subID), whereby each "sufficient" motion bin
    has >= num_sim_per_bin.

    The broad steps are:
    1) Execute `simulate_run` function (which has up to 2 rounds of sampling).
       Sampling = generating simulation runs and binning them in finer bins.
                  For example, in intervals of 10% for mean % censored frames
                  of simulation run.
    2) Execute `select_noncensored_sim_run` function (which first randomly picks
       simulation runs from each sampled finer bin until each fine bin has equal
       number of simulation runs. Then for selection phase, each bin will be
       coarser. For example, in intervals of 25% for mean % censored frames of
       simulation run. For selection coarser bin (0-25%], will randomly select
       from sampled finer bins (0-10%], (10-20%], & simulation runs that are
       (20-25%] in sampled finer bin (20-30%])
       Selection = same as sampling, except binning into coarser bins.
                   Rationale for this is to ensure that coarser bins will
                   have a good mix of mean % censored frames from the
                   finer bins (since each finer bin will be forced to have
                   equal number of samples before selection).

    Results (the select dictionary, the selected windowID dataframe, and
    diagnostic plots) are saved to disk under the simulation result
    directories defined in `pl_config_file`; nothing is returned.

    Args:
        main_config_file (str): Path to the main config JSON file. Note the
            function internally overwrites this with a fixed replication
            path before reading it, so the passed-in value is not actually
            used.
        pl_config_file (str): Path to the pipeline config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
    """
    # Read main config
    main_config_file = os.path.join(PROJECT_DIR, 'replication',
                                    'config', 'main_config.json')
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rep_sel_subID_window_dir = pl_config['path']['rep_sel_subID_window_dir']
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    sample_bin_interval = pl_config['select_subID']['sim'][
        'sample_bin_interval']
    sample_bin_min = pl_config['select_subID']['sim']['sample_bin_min']
    sample_bin_max = pl_config['select_subID']['sim']['sample_bin_max']
    num_sim_per_bin = pl_config['select_subID']['sim']['num_sim_per_bin']
    num_sim_run_per_sim_ses = pl_config['select_subID']['sim'][
        'num_sim_run_per_sim_ses']
    select_bin_min = pl_config['select_subID']['sim']['select_bin_min']
    select_bin_max = pl_config['select_subID']['sim']['select_bin_max']
    select_bin_interval = pl_config['select_subID']['sim'][
        'select_bin_interval']
    rm_run_th = pl_config['select_subID']['sim']['rm_run_th']

    # Assign input files
    windowID_list = os.path.join(rep_sel_subID_window_dir, 'window_list_1.txt')

    # Read input files
    df = pd.read_csv(windowID_list, dtype={'bld': str})

    # Initialize list of sampling bins
    # Note difference between sampling vs selection bins as stated above
    # in `select_subID_sim` function description.
    sample_bin_list = generate_bin_list(sample_bin_min, sample_bin_max,
                                        sample_bin_interval)

    # Extract subID-level df
    df_subID = util.filter_dataset_subID_from_all_subID_df(df, dataset, subID)
    df_subID = df_subID.reset_index(drop=True)

    # Compute number of simulation runs to select from per bin
    target_count_per_bin = num_sim_per_bin * num_sim_run_per_sim_ses

    # Simulate runs to generate `selected_window_idx`,
    # a np array of # windows per run x # simulated runs
    selected_window_idx, sim_run_percent_censored, sim_run_FD_metric, sim_run_DV_metric, sufficient_bins = simulate_run(
        pl_config_file, df_subID, sample_bin_list, target_count_per_bin,
        dataset, subID)

    # Initialize list of selection bins
    # Note difference between sampling vs selection bins as stated above
    # in `select_subID_sim` function description.
    select_bin_list = generate_bin_list(select_bin_min, select_bin_max,
                                        select_bin_interval)

    #############################################################################
    # ===========================================================================
    # PRE-CENSORED BINS
    # ===========================================================================
    #############################################################################
    # Drop select bins that do not have a wide enough spread of
    # % censored frames.
    select_bin_counts = compute_bin_count(sim_run_percent_censored,
                                          select_bin_list)

    # Drop bins with 0 counts
    sufficient_select_bins = filter_zero_count_bin(pl_config_file,
                                                   select_bin_counts,
                                                   sim_run_percent_censored)

    #############################################################################
    # ===========================================================================
    # CENSORED BINS
    # ===========================================================================
    #############################################################################
    # Generate list of sufficient select bins
    # (excluding non-censored bin)
    censored_bin_list = []
    for bin_idx, bin_count in sufficient_select_bins.items():
        bin_min, bin_max = bin_idx.left, bin_idx.right
        if bin_max != 0:
            censored_bin_list.append(f"{bin_min}_{bin_max}")

    # Determine combinations of selection (coarser) bins &
    # number of simulation runs (within each simulation
    # session)
    #
    # Imagine there are 4 selection bins & 2 simulation runs
    # (within each simulation session). This would form a
    # 4x4 grid (because each simulation run has 4 selection
    # bins). Only retain combinations (not permutations).
    # This means that we only take the bottom right
    # triangle (inclusive of diagonal) of the 4x4 grid.
    combo_dict = get_combination(pl_config_file, censored_bin_list)

    # Init dictionary storing all final analysis
    # % censored frames threshold (to discard run) &
    # combinations of bins
    # Hierarchy is:
    #   rm_run_th_value
    #       a list of lists; each nested list comprises
    #       a combination of the bins.
    select_dict = initialize_select_dict(pl_config_file, combo_dict)

    # Select censored simulation runs for each bin for each instance.
    # Only executes for the loosest (to discard run) % censored frames
    # threshold.
    select_dict = select_censored_bin(pl_config_file, combo_dict, select_dict,
                                      sim_run_percent_censored,
                                      selected_window_idx, dataset, subID)

    #############################################################################
    # ===========================================================================
    # NON-CENSORED BINS
    # ===========================================================================
    #############################################################################
    # Select 0% censored frames instances
    for bin_idx, bin_count in sufficient_select_bins.items():
        bin_min, bin_max = bin_idx.left, bin_idx.right
        if bin_max == 0:
            noncensored_bin_min, noncensored_bin_max = bin_min, bin_max
            bin_key = f"{noncensored_bin_min}_{noncensored_bin_max}"
            select_dict, select_noncensored_bin_success = select_noncensored_bin(
                pl_config_file, dataset, subID, bin_key, select_dict,
                selected_window_idx, sim_run_percent_censored)

            # Append 0% censored frames combinations to `combo_dict`
            if select_noncensored_bin_success:
                loosest_rm_run_th_value = np.max(rm_run_th)
                noncensored_bin_list = []
                for _ in range(num_sim_run_per_sim_ses):
                    noncensored_bin_list.append(
                        f"{noncensored_bin_min}_{noncensored_bin_max}")
                    combo_dict[loosest_rm_run_th_value].append(
                        noncensored_bin_list.copy())

    #############################################################################
    # ===========================================================================
    # POST-NON-CENSORED BINS
    # ===========================================================================
    #############################################################################
    # Save `select_dict`
    select_dict_dir = os.path.join(rep_sel_subID_sim_dir, 'result', 'selected',
                                   'select_dict')
    os.makedirs(select_dict_dir, exist_ok=True)
    select_dict_file = os.path.join(select_dict_dir, f"{dataset}_{subID}.pkl")
    with open(select_dict_file, "wb") as f:
        pickle.dump(select_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Drop non-unique bin combinations to avoid repeated analyses.
    combo_dict = drop_non_unique_combination(combo_dict)
    windowID_df = generate_sim_run_df(pl_config_file, df_subID, combo_dict,
                                      select_dict)

    # Re-order output df columns
    first_column = [
        'dataset', 'subID', 'rm_run_th', 'combo', 'sim_run_bin', 'sim_sesID',
        'sim_runID', 'sim_bld', 'windowID', 'window_idx_start',
        'window_idx_end', 'bld', 'sim_run_percent_censored',
        f'sim_run_{FD_metric}', f'sim_run_{DV_metric}'
    ]
    windowID_df = util.reorder_df_columns(windowID_df, first_column)

    # Save windowID_df
    selected_df_dir = os.path.join(rep_sel_subID_sim_dir, 'result', 'selected',
                                   'df')
    filename = f"{dataset}_{subID}.csv"
    util.save_df_to_csv(windowID_df, selected_df_dir, filename)
    print(
        f"Simulation run-level csv for selected {dataset} {subID} saved to {os.path.join(selected_df_dir, filename)}"
    )

    # Plot selected subID mean % censored frames distribution
    selected_dir = os.path.join(rep_sel_subID_sim_dir, 'result', 'selected')
    windowID_df['sim_run_percent_censored_first'] = windowID_df.groupby(
        ['dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
         'sim_runID'])['sim_run_percent_censored'].transform('first')
    windowID_compressed_df = windowID_df.drop_duplicates(subset=[
        'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID', 'sim_runID'
    ])
    plot_distribution_of_sim_run_percent_censored(
        df_subID=windowID_compressed_df,
        selected_window_idx=None,
        output_dir=selected_dir,
        dataset=dataset,
        subID=subID,
        bin_interval=5,
        color='green',
        column_name='sim_run_percent_censored')

    # Plot selected subID heatmap for 2 simulation runs
    # within 1 simulation session set-up
    if num_sim_run_per_sim_ses == 2:
        plot_1ses2run_no_rm_run_th_heatmap(pl_config_file, select_dict,
                                           dataset, subID)


# Wrapper 2
def simulate_run(pl_config_file,
                 df_subID,
                 bin_list,
                 target_count_per_bin,
                 dataset,
                 subID,
                 examples=False):
    """Simulate fMRI runs from non-ground-truth windows, using up to 2 rounds of sampling.

    Overall, this function simulates fMRI runs from non-ground truth
    windows.

    The broad steps are:
    1) 1st round sampling of windows.
        - Weighted sampling performed.
        - Windows from rarer motion bins are weighted with higher
          probability of selection.
        - Rationale is so don't need as many sampling iterations to
          simulate that motion bin.
    2) If subID has bins with insufficient simulated runs, then
       sample more for those "deficient" bins (2nd round sampling).
       - This is done in a "up" or "down" fashion.
       - Example of "up": if subID missing bins (70-80%],
         (80-90%], & (90-100%], then when sampling for (70-80%]
         bin, will sample windows from bins (70-80%], (80-90%],
         & (90-100%]. When sampling for (80-90%], will sample
         windows from bins (80-90%] & (90-100%]. When sampling
         for (90-100%], will only sample from (90-100%] bin.
       - Example of "down": if subID missing bins (0-10%], & (10-20%],
         then when sampling for (10-20%], will sample windows from
         bins (0-10%], & (10-20%]. When sampling for (0-10%], will
         sample windows from only bin (0-10%].
       - Rationale for "up" & "down" sampling is because these deficient
         bins tend to be at the extremes (e.g., (0-10%] or (90-100%]).
         So randomly sampling from non-extreme bins will require many
         iterations for all selected windows to be from the same bin.
         Thus, this allows for more efficient sampling.

    After sampling, simulation runs that are non-unique, contain repeated
    windowIDs, or are fully (100%) censored are filtered out. If not called
    in `examples` mode and the candidate subID ends up with fewer than
    `num_sample_bin2select_subID` sufficient censored bins, `sim_reject_subID`
    is invoked to terminate the algorithm for this subID.

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        df_subID (pd.DataFrame): Window-level dataframe for a single subID.
        bin_list (list): List of bin edges (e.g. from `generate_bin_list`)
            used for sampling-round binning of % censored frames.
        target_count_per_bin (int): Minimum number of simulation run
            instances required per bin to be considered "sufficient".
        dataset (str): Dataset name.
        subID (str): Subject ID.
        examples (bool, optional): If True, skips candidate-plot generation
            and the insufficient-bin rejection check (used when calling this
            function for illustrative/example purposes rather than as part
            of the full selection pipeline). Defaults to False.

    Returns:
        tuple: A 5-tuple of
            (selected_window_idx (np.ndarray): window indices for each
            simulated run, shape (num_window_per_run, num simulated runs);
            sim_run_percent_censored (np.ndarray): mean % censored frames
            per simulation run; sim_run_FD_metric (np.ndarray): mean FD
            metric per simulation run; sim_run_DV_metric (np.ndarray): mean
            DV metric per simulation run; sufficient_bins (pd.Series): bin
            counts for bins meeting `target_count_per_bin`).
    """
    print(f"Start simulating runs for {dataset} {subID}...")

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    sim_iter = pl_config['select_subID']['sim']['sim_iter']
    sample_bin_min = pl_config['select_subID']['sim']['sample_bin_min']
    num_window_per_run = pl_config['select_subID']['sim']['num_window_per_run']
    num_sample_bin2select_subID = pl_config['select_subID']['sim'][
        'num_sample_bin2select_subID']
    sim_deficient_bin_iter = pl_config['select_subID']['sim'][
        'sim_deficient_bin_iter']
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']

    # Weighted sampling of windows
    subID_bin_prob_series = compute_bin_prob(df_subID, bin_list)
    window_prob_arr = compute_window_prob(subID_bin_prob_series)
    window_idx = df_subID.index.to_numpy()
    selected_window_idx = sample_weighted_window(window_idx,
                                                 num_window_per_run, sim_iter,
                                                 window_prob_arr)

    # Filter simulation runs
    selected_window_idx = filter_wrapper(selected_window_idx)

    # Identify bins with insufficient simulated runs
    sim_run_percent_censored = compute_column_mean_of_sim_run_from_window(
        df_subID, selected_window_idx, column_name='window_percent_censored')
    bin_counts = compute_bin_count(sim_run_percent_censored, bin_list)
    deficient_bins = bin_counts[bin_counts < target_count_per_bin]

    # If there are bins with insufficient simulated runs, sample more
    if not deficient_bins.empty:
        print(f"Simulating more runs for {dataset} {subID}...")
        selected_window_idx = sample_more_for_deficient_bins(
            deficient_bins, bin_list, df_subID, sample_bin_min,
            num_window_per_run, sim_deficient_bin_iter, selected_window_idx)
        # Filter simulation runs
        selected_window_idx = filter_wrapper(selected_window_idx)

    # Filter out simulation runs with 100% censored frames
    sim_run_percent_censored = compute_column_mean_of_sim_run_from_window(
        df_subID, selected_window_idx, column_name='window_percent_censored')
    selected_window_idx, sim_run_percent_censored = filter_fully_censored_sim_run(
        selected_window_idx, sim_run_percent_censored)

    # Identify bins with sufficient simulated runs
    bin_counts = compute_bin_count(sim_run_percent_censored, bin_list)
    sufficient_bins = bin_counts[bin_counts >= target_count_per_bin]

    # Compute other simulation run-level summary metrics
    sim_run_FD_metric = compute_column_mean_of_sim_run_from_window(
        df_subID, selected_window_idx, column_name=f'window_mean_{FD_metric}')
    sim_run_DV_metric = compute_column_mean_of_sim_run_from_window(
        df_subID, selected_window_idx, column_name=f'window_mean_{DV_metric}')

    # Plot candidate subID mean % censored frames distribution
    if not examples:
        rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']
        candidate_dir = os.path.join(rep_sel_subID_sim_dir, 'result',
                                     'candidate')
        plot_distribution_of_sim_run_percent_censored(
            df_subID,
            selected_window_idx,
            candidate_dir,
            dataset,
            subID,
            bin_interval=10,
            color='blue',
            column_name='window_percent_censored')

        # If candidate subID has insufficient bins, terminate algorithm
        # Only count censored bins (i.e., exclude non-censored bin)
        censored_sufficient_bins = sufficient_bins.loc[
            ~sufficient_bins.index.astype(str).str.contains("-")]
        sufficient_bin_count = len(censored_sufficient_bins)
        if sufficient_bin_count < num_sample_bin2select_subID:
            sim_reject_subID(pl_config_file, sufficient_bin_count, dataset,
                             subID, df_subID, selected_window_idx)

    return selected_window_idx, sim_run_percent_censored, sim_run_FD_metric, sim_run_DV_metric, sufficient_bins


# Wrapper 3
def filter_wrapper(selected_window_idx):
    """Filter selected window indices to remove non-unique and internally-repeated simulation runs.

    Filters the selected window indices based on two criteria:
    1. Removes non-unique simulation runs.
    2. Removes simulation runs that have repeated windowIDs within simulation run.

    Args:
        selected_window_idx (np.ndarray): The indices of the selected windows
            for simulation runs, shape (num_window_per_run, num simulated
            runs).

    Returns:
        np.ndarray: The filtered window indices after removing non-unique
            simulation runs and removing repeated windowIDs within
            simulation runs.
    """

    # Step 1: Filter out non-unique simulation runs
    selected_window_idx = filter_unique_sim_run(selected_window_idx)

    # Step 2: Filter out simulation runs with repeated windowID (within the simulation run)
    selected_window_idx = filter_sim_run_with_repeated_windowID(
        selected_window_idx)

    return selected_window_idx


# Wrapper 4
def select_censored_bin(pl_config_file, combo_dict, select_dict,
                        sim_run_percent_censored, selected_window_idx, dataset,
                        subID):
    """Select censored simulation run instances for every run/bin combination in `select_dict`.

    For the 1st simulation run, instances are chosen simultaneously across
    all of its bins by weighted random selection (weighting is toward the
    edges of each bin's % censored frames range via
    `compute_weighted_probability_vector`). For the 2nd and subsequent
    simulation runs, bins are processed one at a time (most "problematic"
    bin, i.e. sorted last alphabetically/numerically, first), and instances
    are drawn one-at-a-time so that, for a given instance index, the windows
    selected for this run's bin do not overlap with the windows already
    selected for previous runs' bins at that same instance index (windows
    cannot repeat across simulation runs belonging to the same simulation
    session/instance).

    This selection process is retried up to `num_select_try` times (with a
    different random seed each time) until every run/bin combination in
    `select_dict` is filled. If no attempt fully succeeds, the attempt with
    the fewest unfilled bins is tracked, and its seed is reused on the final
    attempt to give the best achievable result.

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        combo_dict (dict): Dictionary (keyed by rm_run_th value) of lists of
            bin-combination tuples/lists, as produced by `get_combination`.
        select_dict (dict): Nested dictionary of {run_key: {bin_key: None}}
            (as produced by `initialize_select_dict`) to be filled in with
            selected window-index arrays.
        sim_run_percent_censored (np.ndarray): Mean % censored frames for
            each candidate simulation run instance.
        selected_window_idx (np.ndarray): Window indices for each candidate
            simulation run instance, shape (num_window_per_run, num
            candidate instances).
        dataset (str): Dataset name.
        subID (str): Subject ID.

    Returns:
        dict: `select_dict`, updated in place so that each run/bin
            combination maps to a np.ndarray of selected window indices
            (shape (num_window_per_run, num_sim_per_bin)) where selection
            succeeded, or remains None where it did not.
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    num_select_try = pl_config['select_subID']['sim']['num_select_try']
    pre_select_bin_interval = pl_config['select_subID']['sim'][
        'pre_select_bin_interval']
    num_sim_per_bin = pl_config['select_subID']['sim']['num_sim_per_bin']
    rm_run_th = pl_config['select_subID']['sim']['rm_run_th']

    # Init variables
    select_success = False
    previous_num_insufficient_bin = 99999
    best_select_try_idx = None

    for select_try_idx in range(num_select_try):
        print(
            f"Attempt {select_try_idx + 1}/{num_select_try} to select censored simulation run instances..."
        )

        # Init variables
        first_run = True
        remain_instance_window = selected_window_idx
        remain_instance_percent_censored = sim_run_percent_censored
        loosest_rm_run_th = np.max(rm_run_th)

        # Loop through each simulation run ID
        for sim_run_idx, sim_runID in enumerate(list(select_dict.keys())):
            # Define run-level variables
            run_key = f"run-{str((sim_run_idx + 1)).zfill(2)}"
            run_bins = list(select_dict[sim_runID].keys())

            # If select instances for 1st run, weighted selection across each first run's bins simultaneously.
            if first_run:
                for run_bin_idx, run_bin in enumerate(run_bins):
                    # If final try, apply best seed.
                    if (select_try_idx
                            == num_select_try - 1) and (best_select_try_idx):
                        np.random.seed(run_bin_idx + best_select_try_idx)
                    # Else try another seed.
                    else:
                        np.random.seed(run_bin_idx + select_try_idx)

                    # Extract indices of instances within bin range
                    run_bin_min, run_bin_max = map(int, run_bin.split("_"))
                    run_bin_idx = np.where(
                        (remain_instance_percent_censored > run_bin_min)
                        & (remain_instance_percent_censored <= run_bin_max))[0]
                    remain_instance_bin_percent_censored = remain_instance_percent_censored[
                        run_bin_idx]

                    # Compute weighted probability to select each instance
                    p_vector = compute_weighted_probability_vector(
                        run_bin_min, run_bin_max, pre_select_bin_interval,
                        remain_instance_bin_percent_censored)

                    # Weighted selection of instances across bin
                    select_instance_idx = np.random.choice(
                        run_bin_idx,
                        size=num_sim_per_bin,
                        replace=False,
                        p=p_vector)

                    # Map indices of instances back to get windows of instances
                    select_instance_window = remain_instance_window[:,
                                                                    select_instance_idx]

                    # Store instances in dictionary
                    select_dict[run_key][run_bin] = select_instance_window

                    # Update remaining instances for subsequent loops
                    remain_instance_percent_censored = np.delete(
                        remain_instance_percent_censored, select_instance_idx)
                    remain_instance_window = np.delete(remain_instance_window,
                                                       select_instance_idx,
                                                       axis=1)

                # Set false for 2nd run onwards
                first_run = False

            # Else, weighted selection for each run for each bin, 1-instance-at-a-time.
            else:
                # Sort run_bins starting with "problematic" bin with highest probability
                # of overlapping windows with previous runs.
                sorted_run_bins = sorted(run_bins, reverse=True)

                # Loop through each bin
                for sorted_run_bin in sorted_run_bins:
                    # Determine bin limits
                    run_bin_min, run_bin_max = map(int,
                                                   sorted_run_bin.split("_"))

                    # Identify previous runs' bins which are associated with current
                    # run's bin to avoid overlapping windows.
                    previous_instance_window_list = []
                    seen_run_combinations = set()

                    for combo in combo_dict[loosest_rm_run_th]:
                        if combo[sim_run_idx] == sorted_run_bin:
                            for previous_run_idx, previous_run_bin in enumerate(
                                    combo[:sim_run_idx]):
                                previous_run_key = f"run-{str((previous_run_idx + 1)).zfill(2)}"
                                run_combination = (previous_run_key,
                                                   previous_run_bin)
                                if run_combination not in seen_run_combinations:
                                    seen_run_combinations.add(run_combination)
                                    previous_run_bin_arr = select_dict[
                                        previous_run_key][previous_run_bin]
                                    if isinstance(previous_run_bin_arr,
                                                  np.ndarray):
                                        previous_instance_window_list.append(
                                            previous_run_bin_arr)

                    # Concatenate runs' bins' (to-avoid-conflict-with) windows
                    previous_instance_window = np.vstack(
                        previous_instance_window_list)

                    # Loop through each instance
                    select_instance_window_list = []
                    for instance_idx in range(num_sim_per_bin):
                        # Set random seed
                        np.random.seed(instance_idx)

                        # Check sufficient instances to select
                        if np.shape(
                                remain_instance_window)[1] < num_sim_per_bin:
                            raise ValueError(
                                f"{dataset} {subID} {sorted_run_bin} insufficient "
                                f"({np.shape(remain_instance_window)[1]}/{num_sim_per_bin}) simulation run "
                                f"instances for {sim_run_idx + 1}-th simulation run"
                            )

                        # Extract corresponding instances from previous runs' bins to avoid overlapping windows.
                        previous_instance_window_avoid = previous_instance_window[:,
                                                                                  instance_idx]
                        overlap_window_bool = np.isin(
                            remain_instance_window,
                            previous_instance_window_avoid)
                        no_overlap_window_col = np.all(~overlap_window_bool,
                                                       axis=0)
                        no_overlap_window_idx = np.where(
                            no_overlap_window_col)[0]

                        # Extract indices of instances within bin range
                        run_bin_idx = no_overlap_window_idx[np.where(
                            (remain_instance_percent_censored[
                                no_overlap_window_idx] > run_bin_min)
                            & (remain_instance_percent_censored[
                                no_overlap_window_idx] <= run_bin_max))[0]]

                        if run_bin_idx.size == 0:
                            print(
                                f"WARNING: {dataset} {subID} {sim_runID} {sorted_run_bin} index {instance_idx} "
                                f"instance had {np.shape(remain_instance_window)[1]} candidate windows to select "
                                f"from, but these windows have >= 1 overlapping window with previous "
                                f"simulation runs, no windows left to select. Proceeding with next bin..."
                            )
                            break

                        remain_instance_bin_percent_censored = remain_instance_percent_censored[
                            run_bin_idx]

                        # Compute weighted probability to select instance
                        p_vector = compute_weighted_probability_vector(
                            run_bin_min, run_bin_max, pre_select_bin_interval,
                            remain_instance_bin_percent_censored)

                        # Weighted selection of instance across bin
                        select_instance_idx = np.random.choice(run_bin_idx,
                                                               size=1,
                                                               replace=False,
                                                               p=p_vector)

                        # Map indices of instances back to get windows of instances
                        select_instance_window = remain_instance_window[:,
                                                                        select_instance_idx]

                        # Append instance in list
                        select_instance_window_list.append(
                            select_instance_window)

                        # Update remaining instances for subsequent loops
                        remain_instance_percent_censored = np.delete(
                            remain_instance_percent_censored,
                            select_instance_idx)
                        remain_instance_window = np.delete(
                            remain_instance_window,
                            select_instance_idx,
                            axis=1)

                    # Concatenate all instances within run
                    if len(select_instance_window_list) == num_sim_per_bin:
                        select_dict[run_key][sorted_run_bin] = np.hstack(
                            select_instance_window_list)

        # Successful if all runs' bins are filled.
        select_success = all(
            isinstance(arr, np.ndarray) for run_key in select_dict
            for arr in select_dict[run_key].values())

        # If successful, stop selecting prematurely.
        if select_success:
            print(f"Succeeded to select censored simulation run instances on "
                  f"Attempt {select_try_idx + 1}/{num_select_try}.")
            break
        else:
            num_insufficient_bin = sum(
                not isinstance(arr, np.ndarray) for run_key in select_dict
                for arr in select_dict[run_key].values())
            if num_insufficient_bin < previous_num_insufficient_bin:
                previous_num_insufficient_bin = num_insufficient_bin
                best_select_try_idx = select_try_idx

    return select_dict


# Wrapper 5
def select_noncensored_bin(pl_config_file, dataset, subID, bin_key,
                           select_dict, selected_window_idx,
                           sim_run_percent_censored):
    """Select 0%-censored ("non-censored") simulation run instances for every simulation run.

    For the 1st simulation run, `num_sim_per_bin` instances with 0% mean
    censored frames are randomly selected (vectorized). For the 2nd and
    subsequent simulation runs, instances are drawn one-at-a-time so that,
    for a given instance index, the windows selected for this run do not
    overlap with the windows already selected for previous runs at that same
    instance index (windows cannot repeat across simulation runs belonging
    to the same simulation session/instance). Unlike `select_censored_bin`,
    this function does not retry with different seeds: if at any point there
    are insufficient 0%-censored instances (fewer than `num_sim_per_bin`)
    left to choose from, or none of the remaining candidates avoid
    overlapping with previous runs, selection is aborted immediately and the
    unmodified `select_dict` is returned together with a False success flag.

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        bin_key (str): Key (e.g. "0_0") identifying the non-censored bin
            within `select_dict`.
        select_dict (dict): Nested dictionary of {run_key: {bin_key: ...}}
            to be filled in with selected window-index arrays for the
            non-censored bin.
        selected_window_idx (np.ndarray): Window indices for each candidate
            simulation run instance, shape (num_window_per_run, num
            candidate instances).
        sim_run_percent_censored (np.ndarray): Mean % censored frames for
            each candidate simulation run instance.

    Returns:
        tuple: A 2-tuple of (select_dict (dict): a deep copy of the input
            `select_dict` with the non-censored bin filled in for every
            simulation run if successful, otherwise the original
            `select_dict` unchanged; success (bool): True if the
            non-censored bin was successfully filled for every simulation
            run, False otherwise).
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    num_sim_per_bin = pl_config['select_subID']['sim']['num_sim_per_bin']
    num_sim_run_per_sim_ses = pl_config['select_subID']['sim'][
        'num_sim_run_per_sim_ses']

    # Copy `select_dict`
    out_select_dict = copy.deepcopy(select_dict)

    # Loop through each simulation run
    for sim_run_idx in range(num_sim_run_per_sim_ses):
        # Define run_key
        run_key = f"run-{str((sim_run_idx + 1)).zfill(2)}"

        # If 1st simulation run, vectorized sampling of instances
        # of 1st simulation run
        if sim_run_idx == 0:
            # Identify indices whereby mean % censored frames = 0
            # (i.e., non-censored instances)
            instance_bin_idx = np.where(sim_run_percent_censored == 0)[0]

            # Check sufficient instances to sample
            if len(instance_bin_idx) < num_sim_per_bin:
                print(
                    f"WARNING: {dataset} {subID} 0% censored bin insufficient "
                    f"({len(instance_bin_idx)}/{num_sim_per_bin}) simulation run "
                    f"instances for {sim_run_idx + 1}-th simulation run. "
                    f"Will skip selecting 0% censored bin...")
                return select_dict, False
            # Randomly select 0% mean censored frames instances
            select_instance_bin_idx = np.random.choice(instance_bin_idx,
                                                       num_sim_per_bin,
                                                       replace=False)

            # Map indices of instances back to get windows of instances
            select_instance_bin_window = selected_window_idx[:,
                                                             select_instance_bin_idx]

            # Append selected instance windows to dictionary
            out_select_dict[run_key][bin_key] = select_instance_bin_window

            # Discard used instances
            remain_instance_bin_idx = np.setdiff1d(instance_bin_idx,
                                                   select_instance_bin_idx)
            remain_instance_bin_window = selected_window_idx[:,
                                                             remain_instance_bin_idx]

            # Check no missing instances
            num_bin_start_instance = len(instance_bin_idx)
            num_remain_instance = np.shape(remain_instance_bin_window)[1]
            num_select_instance = np.shape(select_instance_bin_window)[1]
            assert num_bin_start_instance == num_remain_instance + num_select_instance, (
                f"ERROR: when selecting instances for 0% censored bin for 1st simulation run, "
                f"num_bin_start_instance != num_remain_instance + num_select_instance. Obtained "
                f"{num_bin_start_instance} == {num_remain_instance} + {num_select_instance}."
            )

        # If not 1st simulation run, need to loop through each
        # instance of N-th simulation run to make sure no
        # overlapping windows with all previous simulation runs
        # for that n-th instance.
        #
        # For example, the 5th instance for
        # 1st simulation run & 5th instance for 2nd simulation run
        # (this means both 1st & 2nd simulation runs under the same
        # simulation session) cannot have overlapping windows.
        else:
            # Concatenate selected instances' windows vertically, within-bin,
            # across all previous simulation runs.
            previous_sim_run_window_list = []
            for previous_sim_run_idx in range(sim_run_idx):
                # Append previous_run_key
                previous_run_key = f"run-{str((previous_sim_run_idx + 1)).zfill(2)}"
                previous_sim_run_window_list.append(
                    out_select_dict[previous_run_key][bin_key])
            previous_sim_run_window = np.vstack(previous_sim_run_window_list)

            # Loop through each instance to only randomly sample
            # candidate instances with non-overlapping windows
            for instance_idx in range(num_sim_per_bin):
                # Check sufficient instances to sample
                if np.shape(remain_instance_bin_window)[1] < num_sim_per_bin:
                    print(
                        f"{dataset} {subID} 0% censored bin insufficient "
                        f"({np.shape(remain_instance_bin_window)[1]}/{num_sim_per_bin}) simulation run "
                        f"instances for {sim_run_idx + 1}-th simulation run.")
                    return select_dict, False

                # Check for instances' indices without overlapping windows
                previous_sim_run_instance_window = previous_sim_run_window[:,
                                                                           instance_idx]
                overlap_window_bool = np.isin(
                    remain_instance_bin_window,
                    previous_sim_run_instance_window)
                no_overlap_window_col = np.all(~overlap_window_bool, axis=0)
                no_overlap_window_idx = np.where(no_overlap_window_col)[0]

                # Check if no_overlap_window_idx is empty
                if no_overlap_window_idx.size == 0:
                    print(
                        f"Due to non-overlapping windows constraint, "
                        f"{dataset} {subID} 0% censored bin insufficient "
                        f"({np.shape(remain_instance_bin_window)[1]}/{num_sim_per_bin}) simulation run "
                        f"instances for {sim_run_idx + 1}-th simulation run.")
                    return select_dict, False
                else:
                    # Randomly select 0% mean censored frames instance
                    select_instance_bin_idx = np.random.choice(
                        no_overlap_window_idx, 1, replace=False)

                    # Map indices of instances back to get remaining windows of instances
                    select_instance_bin_window = remain_instance_bin_window[:,
                                                                            select_instance_bin_idx]

                    # Double-check no common window IDs between n-th instance of N-th simulation run &
                    # n-th instance of all previous (N-1)-th, (N-2)-th, ..., 2nd, 1st, simulation run
                    common_window = np.intersect1d(
                        select_instance_bin_window,
                        previous_sim_run_instance_window)
                    assert len(common_window) == 0, (
                        f"ERROR: select_instance_bin_window & previous_sim_run_instance_window "
                        f"contain {len(common_window)} common window IDs.")

                    # Append selected instance windows to dictionary
                    if bin_key not in out_select_dict[run_key]:
                        out_select_dict[run_key][
                            bin_key] = select_instance_bin_window
                    else:
                        out_select_dict[run_key][bin_key] = np.hstack(
                            (out_select_dict[run_key][bin_key],
                             select_instance_bin_window))
                    # Discard used instances
                    remain_instance_bin_window = np.delete(
                        remain_instance_bin_window,
                        select_instance_bin_idx,
                        axis=1)

    return out_select_dict, True


#############################################################################
# ===========================================================================
# UTIL FUNCTIONS
# ===========================================================================
#############################################################################
def generate_bin_list(bin_min, bin_max, bin_interval):
    """Generate a list of bin edges based on the provided minimum, maximum, and interval values.

    Example:
        bin_min = -10
        bin_max = 100
        bin_interval = 10
        bin_list = generate_bin_list(bin_min, bin_max, bin_interval)
        bin_list value is [-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    Args:
        bin_min (int): Minimum bin edge value.
        bin_max (int): Maximum bin edge value.
        bin_interval (int): Interval (width) between successive bin edges.

    Returns:
        list: List of int bin edges from `bin_min` to `bin_max` (inclusive),
            spaced by `bin_interval`.
    """
    return list(range(bin_min, bin_max + bin_interval, bin_interval))


def compute_bin_prob(df_subID, bin_list):
    """Compute the (inverse-frequency-weighted) selection probability for each window's bin.

    Weighting is performed. For example, if a window is in the
    (0, 10]% bin, and (0, 10]% bin windows are very rare for a
    certain subID, then there will be a higher probability of such
    a rare window being selected.

    Args:
        df_subID (pd.DataFrame): Window-level dataframe for a single subID;
            must contain a 'window_percent_censored' column. A 'bin' column
            is added to this dataframe as a side effect.
        bin_list (list): List of bin edges used to categorize
            'window_percent_censored' via `pd.cut`.

    Returns:
        pd.Series: A Series (aligned with `df_subID`'s index) containing the
            computed, normalized bin probability for each window.
    """
    df_subID['bin'] = pd.cut(df_subID['window_percent_censored'],
                             bins=bin_list,
                             right=True)
    bin_counts = df_subID['bin'].value_counts(normalize=True)

    # If bin present, take inverted probability.
    # Else, inverted probability = 0 (to avoid error where 1 divided by 0 = inf)
    bin_inverted_probs = (1 / bin_counts).where(bin_counts > 0, 0)
    bin_inverted_probs /= bin_inverted_probs.sum()
    bin_inverted_probs_dict = bin_inverted_probs.to_dict()
    bin_prob_series = df_subID['bin'].map(bin_inverted_probs_dict).astype(
        'float64')

    return bin_prob_series


def compute_window_prob(bin_prob_series):
    """Compute normalized per-window selection probabilities from a series of bin probabilities.

    This function takes a pandas Series `bin_prob_series`, which contains probability values (e.g.,
    from a 'bin_prob' column), normalizes these values so that they sum to 1, and then returns the
    resulting probabilities as a numpy array.

    Args:
        bin_prob_series (pd.Series): Series of (unnormalized) per-window bin
            probabilities, e.g. as returned by `compute_bin_prob`.

    Returns:
        np.ndarray: A numpy array of normalized window probabilities, where the sum of the values is 1.
    """

    # Normalize the bin probabilities so that they sum to 1
    window_prob_series = bin_prob_series / bin_prob_series.sum()

    # Convert the normalized probabilities to a numpy array
    window_prob_arr = window_prob_series.values

    return window_prob_arr


def sample_weighted_window(indices, num_window_per_run, sim_iter,
                           window_prob_arr):
    """Randomly sample windows (with replacement) based on weighted probabilities.

    Selects specified number of windows (with replacement) based on
    weighted probabilities.

    Args:
        indices (np.ndarray): Array of candidate window indices to sample
            from.
        num_window_per_run (int): Number of windows to sample per simulated
            run.
        sim_iter (int): Number of simulated runs to generate (i.e. number of
            sampling iterations).
        window_prob_arr (np.ndarray): Array of selection probabilities
            (aligned with `indices`) used for the weighted sampling.

    Returns:
        np.ndarray: A 2D numpy array of shape `(num_window_per_run, sim_iter)`,
                   where each column represents one simulated run
                   ,and each row represents a windowID.
    """
    np.random.seed(0)
    return np.random.choice(indices,
                            size=(num_window_per_run, sim_iter),
                            p=window_prob_arr)


def sort_windowID_column(column):
    """Sort the values of a single simulation run's window-index column and return them as a tuple.

    Used with `np.apply_along_axis` in `filter_unique_sim_run` so that two
    simulation runs containing the same set of windowIDs (regardless of
    sampled order) hash/compare as identical.

    Args:
        column (np.ndarray or array-like): 1D array of window indices for a
            single simulation run.

    Returns:
        tuple: The sorted values of `column` as a tuple.
    """
    return tuple(sorted(column))


def filter_unique_sim_run(window_idx):
    """Filter out non-unique simulation runs.

    Unique simulation run = at least 1 window ID different from other simulation run.

    For example, Simulation Run A contains win-01, win-02, win-03, and
    Simulation Run B contains win-01, win-02, win-03. Then Simulation Run B is not unique,
    and is removed.

    Another example is Simulation Run C contains win-01, win-03, win-02 (sampled in this order), Then
    Simulation Run C is still not unique because it contains the same windowIDs as Simulation Run A.

    Args:
        window_idx (np.ndarray): 2D array of window indices, shape
            (num_window_per_run, num simulated runs), where each column is
            one simulation run.

    Returns:
        np.ndarray: A filtered version of `window_idx` with unique simulation runs
                   (columns with identical sets of window IDs are removed).
    """
    column_sets = np.apply_along_axis(sort_windowID_column,
                                      axis=0,
                                      arr=window_idx)
    _, unique_idx = np.unique(column_sets, axis=1, return_index=True)

    return window_idx[:, unique_idx]


def filter_sim_run_with_repeated_windowID(window_idx):
    """Filter out simulation runs that contain repeated windowIDs within the same run.

    Args:
        window_idx (np.ndarray): 2D array of window indices, shape
            (num_window_per_run, num simulated runs), where each column is
            one simulation run.

    Returns:
        np.ndarray: A filtered version of `window_idx` with only columns
                   that have unique window IDs within each simulation run.
    """
    # Convert window IDs to sets by sorting each column and finding unique row counts
    sorted_window_idx = np.sort(window_idx,
                                axis=0)  # Sort values within each column
    is_unique = np.all(np.diff(sorted_window_idx, axis=0) != 0,
                       axis=0)  # Check for duplicate IDs in columns

    # Filter out columns with duplicates
    filtered_window_idx = window_idx[:, is_unique]

    return filtered_window_idx


def filter_fully_censored_sim_run(window_idx, sim_run_percent_censored):
    """Filter out fully (100%) censored simulation runs (based on strictest censoring).

    Rationale is because if use strictest censoring on this simulation
    run in downstream TMS analysis, will encounter error because there is
    essentially no rsfMRI timeseries frames to work with.

    Assumption is that the same censoring threshold used for censoring
    will be used as censoring threshold for downstream TMS analysis.

    Args:
        window_idx (np.ndarray): 2D array of window indices, shape
            (num_window_per_run, num simulated runs), where each column is
            one simulation run.
        sim_run_percent_censored (np.ndarray): Mean % censored frames for
            each simulation run in `window_idx` (must have the same number
            of entries as `window_idx` has columns).

    Returns:
        tuple: A 2-tuple of (window_idx (np.ndarray): a filtered version of
            `window_idx` with simulation runs that are less than 100%
            censored; sim_run_percent_censored (np.ndarray): the
            correspondingly filtered % censored frames array).
    """

    # Check identical dimension
    assert np.shape(window_idx)[1] == len(sim_run_percent_censored), \
        f"ERROR: `window_idx` has {np.shape(window_idx)[1]} simulation runs, \
            but `sim_run_percent_censored has unexpectedly different \
                {len(sim_run_percent_censored)} simulation runs."

    # Identify indices whereby 100% censored frames
    noncensored_idx = np.where(sim_run_percent_censored == 100)[0]

    # Drop simulation runs
    window_idx = np.delete(window_idx, noncensored_idx, axis=1)
    sim_run_percent_censored = np.delete(sim_run_percent_censored,
                                         noncensored_idx)

    return window_idx, sim_run_percent_censored


def identify_sampling_direction(bin_range, bin_list):
    """Determine the 2nd-round sampling direction ("up" or "down") for a deficient bin.

    Args:
        bin_range (pd.Interval): The deficient bin (e.g. from a
            `pd.cut`/`value_counts` index) whose `.left` edge is looked up
            in `bin_list`.
        bin_list (list): List of bin edges, e.g. from `generate_bin_list`,
            containing the value 50 as one of its edges.

    Returns:
        str: "up" if the bin is above the midpoint (bin containing 50), "down" otherwise.
    """
    bin_idx = bin_list.index(bin_range.left)
    return "up" if bin_idx > bin_list.index(50) else "down"


def sample_more_from_bins(bin_range, bin_list, direction, df_subID,
                          sample_bin_min, num_window_per_run,
                          sim_deficient_bin_iter):
    """Perform 2nd-round sampling from bins in a specified direction ("up" or "down") for a deficient bin.

    Depending on `direction`, all bins from `bin_range` up to the highest bin
    ("up") or all bins from the lowest bin up to `bin_range` ("down") are
    treated as relevant, weighted-sampling is re-run restricted to windows
    from those relevant bins only. See `simulate_run`'s docstring for the
    rationale behind "up"/"down" 2nd-round sampling.

    Args:
        bin_range (pd.Interval): The deficient bin being sampled for.
        bin_list (list): List of bin edges, e.g. from `generate_bin_list`.
        direction (str): "up" or "down", as returned by
            `identify_sampling_direction`.
        df_subID (pd.DataFrame): Window-level dataframe for a single subID.
        sample_bin_min (int): Minimum sampling bin edge. Currently unused
            within this function body.
        num_window_per_run (int): Number of windows to sample per simulated
            run.
        sim_deficient_bin_iter (int): Number of additional simulated runs
            (sampling iterations) to generate for this deficient bin.

    Returns:
        np.ndarray or None: Array of new sampled indices if sampling is
            successful; None if no relevant windows are available for
            sampling.
    """
    bin_idx = bin_list.index(bin_range.left)

    # Determine sampling direction
    if direction == "up":
        relevant_bins = bin_list[bin_idx:]
        relevant_bins = np.array(relevant_bins) - 1e-6
        relevant_bins.tolist()
    else:
        relevant_bins = bin_list[:bin_idx + 1]
        relevant_bins = np.array(relevant_bins) + 1e-6
        relevant_bins.tolist()

    # Filter the original DataFrame for relevant bins
    relevant_windows = df_subID[pd.cut(df_subID["window_percent_censored"],
                                       bins=bin_list,
                                       right=True).isin(relevant_bins)]

    # Extract window idx for sample_weighted_window
    relevant_window_idx = relevant_windows.index.to_numpy()
    # Perform weighted sampling based on `window_prob_arr`
    relevant_bin_prob_arr = compute_bin_prob(relevant_windows, bin_list)
    relevant_window_prob_arr = compute_window_prob(relevant_bin_prob_arr)
    if relevant_window_idx.size > 0:
        new_samples = sample_weighted_window(
            indices=relevant_window_idx,
            num_window_per_run=num_window_per_run,
            sim_iter=sim_deficient_bin_iter,
            window_prob_arr=relevant_window_prob_arr,
        )
        return new_samples
    else:
        return None


def compute_column_mean_of_sim_run_from_window(df, selected_window_idx,
                                               column_name):
    """Compute the per-simulation-run mean of a specified window-level column.

    Args:
        df (pd.DataFrame): Window-level dataframe containing `column_name`.
        selected_window_idx (np.ndarray): 2D array of window indices into
            `df`, shape (num_window_per_run, num simulated runs).
        column_name (str): Name of the column in `df` to average.

    Returns:
        np.ndarray: The mean value of `column_name`, computed per simulation
            run (i.e. per column of `selected_window_idx`).
    """
    # Select the specified column from the DataFrame and use the selected indices
    selected_values = df[column_name].values[selected_window_idx]

    # Compute and return the mean
    return np.mean(selected_values, axis=0)


def plot_distribution_of_sim_run_percent_censored(
        df_subID,
        selected_window_idx,
        output_dir,
        dataset,
        subID,
        bin_interval=10,
        color='blue',
        column_name='window_percent_censored'):
    """Plot and save a histogram of mean % censored frames across simulation runs for a subID.

    For each subID, plots distribution of mean % censored frames
    (within sim_runID) and saves as image.

    Args:
        df_subID (pd.DataFrame): Window-level dataframe for a single subID.
        selected_window_idx (np.ndarray or None): 2D array of window indices
            into `df_subID`, used to compute per-run means via
            `compute_column_mean_of_sim_run_from_window`. If None,
            `df_subID` is expected to already be at the simulation-run level
            (one row per run) and is used directly.
        output_dir (str): Directory to save the output plot image to.
        dataset (str): Dataset name (used in the plot title/filename).
        subID (str): Subject ID (used in the plot title/filename).
        bin_interval (int, optional): Histogram bin width. Defaults to 10.
        color (str, optional): Histogram bar color. Defaults to 'blue'.
        column_name (str, optional): Name of the window-level (or
            run-level, if `selected_window_idx` is None) column to plot the
            distribution of. Defaults to 'window_percent_censored'.
    """
    # Compute sim_run_percent_censored
    sim_run_percent_censored = compute_column_mean_of_sim_run_from_window(
        df_subID, selected_window_idx, column_name)

    # Define the bins
    bin_edges = np.arange(-10, 101, bin_interval)

    # Plot the distribution using seaborn
    plt.figure(figsize=(10, 6))
    sns.histplot(
        sim_run_percent_censored,
        bins=bin_edges,
        kde=False,
        edgecolor='black',
        color=color,  # Use the color argument
        stat='count')  # Use 'count' for absolute frequency

    # Customize the plot
    plt.title(f"{dataset} {subID}, Mean % Censored Frames Distribution",
              fontsize=16)
    plt.xlabel('Mean % Censored Frames (within candidate simulated run)',
               fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(np.arange(0, 101, 10))

    # Define the output file path
    output_file_path = os.path.join(output_dir, f"{dataset}_{subID}.png")

    # Ensure the directory exists
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

    # Save the plot with dpi=300
    plt.savefig(output_file_path, dpi=300)
    plt.close()
    print(f"Plot saved to {output_file_path}")


def save_rejected_subID_summary_metric(
        dataset, subID, sim_runID_min_percent_censored,
        sim_runID_max_percent_censored, sim_runID_min_FD_metric,
        sim_runID_max_FD_metric, sim_runID_min_DV_metric,
        sim_runID_max_DV_metric, rep_sel_subID_sim_dir):
    """Save summary metrics (min, max) for percent-censored, FD, and DV metrics to a CSV file.

    Args:
        dataset (str): Dataset name.
        subID (str): Subject ID.
        sim_runID_min_percent_censored (float): Minimum mean % censored
            frames across simulation runs.
        sim_runID_max_percent_censored (float): Maximum mean % censored
            frames across simulation runs.
        sim_runID_min_FD_metric (float): Minimum mean FD metric across
            simulation runs.
        sim_runID_max_FD_metric (float): Maximum mean FD metric across
            simulation runs.
        sim_runID_min_DV_metric (float): Minimum mean DV metric across
            simulation runs.
        sim_runID_max_DV_metric (float): Maximum mean DV metric across
            simulation runs.
        rep_sel_subID_sim_dir (str): Base directory under which the summary
            CSV file (result/rejected/df/{dataset}_{subID}.csv) is saved.
    """
    # Create a dictionary to hold the summary metrics
    rejected_subID_data = {
        'dataset': [dataset],
        'subID': [subID],
        'sim_run_percent_censored_min': [sim_runID_min_percent_censored],
        'sim_run_percent_censored_max': [sim_runID_max_percent_censored],
        'sim_run_mean_FD_min': [sim_runID_min_FD_metric],
        'sim_run_mean_FD_max': [sim_runID_max_FD_metric],
        'sim_run_mean_DVARS_min': [sim_runID_min_DV_metric],
        'sim_run_mean_DVARS_max': [sim_runID_max_DV_metric]
    }

    # Convert the dictionary to a DataFrame
    df_to_save = pd.DataFrame(rejected_subID_data)

    # Define the file path to save the CSV
    file_path = os.path.join(rep_sel_subID_sim_dir, 'result', 'rejected', 'df',
                             f"{dataset}_{subID}.csv")

    # Ensure the directory exists
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    # Save the DataFrame to a CSV file
    df_to_save.to_csv(file_path, index=False)

    print(
        f"Summary metrics for rejected {dataset} {subID} saved to {file_path}")


def compute_bin_count(data, bin_list):
    """Compute the bin counts for the given data, sorted by the bin edges.

    Args:
        data (array-like): Values (e.g. % censored frames per simulation
            run) to bin via `pd.cut`.
        bin_list (list): List of bin edges to bin `data` into.

    Returns:
        pd.Series: A sorted Series of bin counts, indexed by bin interval.
    """
    # Perform binning using pd.cut
    bin_counts = pd.cut(data, bins=bin_list, right=True).value_counts()

    # Sort the bin counts by bin edges
    bin_counts = bin_counts.sort_index()

    return bin_counts


def sample_more_for_deficient_bins(deficient_bins, bin_list, df_subID,
                                   sample_bin_min, num_window_per_run,
                                   sim_deficient_bin_iter,
                                   selected_window_idx):
    """Wrapper for handling deficient bins by sampling additional simulation runs for them.

    This is a wrapper function for handling deficient bins.
    The broad steps are:
    1) identify 2nd round sampling direction.
    2) sample in that direction.
    3) append more simulated run columns to Numpy array

    Args:
        deficient_bins (pd.Series): Bin counts (indexed by bin interval) for
            bins below `target_count_per_bin`, as identified in
            `simulate_run`.
        bin_list (list): List of bin edges, e.g. from `generate_bin_list`.
        df_subID (pd.DataFrame): Window-level dataframe for a single subID.
        sample_bin_min (int): Minimum sampling bin edge (passed through to
            `sample_more_from_bins`).
        num_window_per_run (int): Number of windows to sample per simulated
            run.
        sim_deficient_bin_iter (int): Number of additional simulated runs
            (sampling iterations) to generate per deficient bin.
        selected_window_idx (np.ndarray): Existing 2D array of window
            indices, shape (num_window_per_run, num simulated runs), to
            append newly-sampled runs to.

    Returns:
        np.ndarray: The updated array of selected window indices after handling deficient bins.
    """
    new_samples = []

    # Loop through each deficient bin
    for bin_range, _ in deficient_bins.items():
        # Identify the sampling direction for the deficient bin
        direction = identify_sampling_direction(bin_range, bin_list)

        # Sample additional windows from the deficient bin
        additional_samples = sample_more_from_bins(bin_range, bin_list,
                                                   direction, df_subID,
                                                   sample_bin_min,
                                                   num_window_per_run,
                                                   sim_deficient_bin_iter)
        # Append the new samples to the list
        if additional_samples is not None:
            new_samples.append(additional_samples)

    # If there are any new samples, concatenate them with the selected window indices
    if new_samples:
        new_samples = np.concatenate(new_samples, axis=1)
        selected_window_idx = np.concatenate(
            [selected_window_idx, new_samples], axis=1)

    return selected_window_idx


def compute_range_across_sim_run(sim_run_array):
    """Compute the min and max values across all elements in a numpy array.

    Args:
        sim_run_array (np.ndarray): Array of per-simulation-run values,
            expected to be 1D (shape (N,)) or 2D with a single row (shape
            (1, N)).

    Returns:
        tuple: A 2-tuple of (min_value, max_value) across all elements of
            `sim_run_array`.
    """
    # Check the shape of the array
    if sim_run_array.ndim == 2 and sim_run_array.shape[0] != 1:
        raise ValueError(
            f"Expected array with one row (shape (1, N)) or one-dimensional array (shape (N,)), "
            f"but got shape {sim_run_array.shape}.")

    # Flatten the array in case it's 2D
    flattened_array = sim_run_array.flatten()

    # Compute min and max
    min_value = np.min(flattened_array)
    max_value = np.max(flattened_array)

    return min_value, max_value


def sim_reject_subID(pl_config_file, sufficient_bin_count, dataset, subID,
                     df_subID, selected_window_idx):
    """Reject a candidate subID by saving its summary metrics and terminating the algorithm.

    Checks if number of sufficient bins meets required threshold.
    If not, function terminates algorithm.

    Computes and saves the min/max mean % censored frames, FD metric, and DV
    metric across the candidate subID's simulation runs (via
    `save_rejected_subID_summary_metric`), then calls `sys.exit(0)` to stop
    the program for this subID.

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        sufficient_bin_count (int): Number of sufficient bins the candidate
            subID achieved (used only for the printed message).
        dataset (str): Dataset name.
        subID (str): Subject ID.
        df_subID (pd.DataFrame): Window-level dataframe for a single subID.
        selected_window_idx (np.ndarray): 2D array of window indices for the
            candidate subID's simulation runs.
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    num_select_bin2select_subID = pl_config['select_subID']['sim'][
        'num_select_bin2select_subID']
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']

    print(
        f"{dataset} {subID} only has {sufficient_bin_count}/{num_select_bin2select_subID} "
        f"required bins. Terminating algorithm.")

    # Compute mean stats (within each simulated run)
    sim_runID_mean_percent_censored = compute_column_mean_of_sim_run_from_window(
        df_subID, selected_window_idx, column_name='window_percent_censored')
    sim_runID_mean_FD_metric = compute_column_mean_of_sim_run_from_window(
        df_subID, selected_window_idx, column_name=f'window_mean_{FD_metric}')
    sim_runID_mean_DV_metric = compute_column_mean_of_sim_run_from_window(
        df_subID, selected_window_idx, column_name=f'window_mean_{DV_metric}')

    # Extract min & max stats (across simulated runs)
    sim_runID_min_percent_censored, sim_runID_max_percent_censored = compute_range_across_sim_run(
        sim_runID_mean_percent_censored)
    sim_runID_min_FD_metric, sim_runID_max_FD_metric = compute_range_across_sim_run(
        sim_runID_mean_FD_metric)
    sim_runID_min_DV_metric, sim_runID_max_DV_metric = compute_range_across_sim_run(
        sim_runID_mean_DV_metric)

    # Compute and save the rejected subID summary metric
    save_rejected_subID_summary_metric(
        dataset, subID, sim_runID_min_percent_censored,
        sim_runID_max_percent_censored, sim_runID_min_FD_metric,
        sim_runID_max_FD_metric, sim_runID_min_DV_metric,
        sim_runID_max_DV_metric, rep_sel_subID_sim_dir)

    # Terminate the program
    sys.exit(0)


def generate_sim_run_df(pl_config_file, df_subID, combo_dict, select_dict):
    """Expand the selected instances in `select_dict` into a window-level dataframe.

    Generates a dataframe whereby 1 row = 1 windowID.
    Multiple windowIDs belong to a single instance.

    This dataframe is used to run downstream analyses steps
    such as generating the simulation run rsfMRI NIFTI by
    concatenating windows for a single instance.

    For each % censored frames threshold (`rm_run_th_value`) and each unique
    bin combination under it, loops through each simulation session
    instance (`sim_sesID`) and each simulation run/bin within the
    combination, looks up the corresponding selected windows in
    `select_dict`, and appends the corresponding rows of `df_subID` (in
    window order) to the output dataframe, tagging them with identifying
    columns (rm_run_th, combo, sim_sesID, sim_runID, sim_bld, sim_run_bin,
    and the run-level % censored frames / FD / DV metrics).

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        df_subID (pd.DataFrame): Window-level dataframe for a single subID;
            must contain a 'windowID' column of the form 'win-NN'.
        combo_dict (dict): Dictionary (keyed by rm_run_th value) of lists of
            unique bin-combination tuples/lists, as produced by
            `drop_non_unique_combination`.
        select_dict (dict): Nested dictionary of {run_key: {bin_key:
            selected window-index array}}, as filled in by
            `select_censored_bin`/`select_noncensored_bin`.

    Returns:
        pd.DataFrame: Window-level dataframe (1 row per windowID per
            instance) with the identifying/summary columns described above
            appended.
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    select_bin_min = pl_config['select_subID']['sim']['select_bin_min']
    select_bin_max = pl_config['select_subID']['sim']['select_bin_max']
    select_bin_interval = pl_config['select_subID']['sim'][
        'select_bin_interval']
    num_sim_per_bin = pl_config['select_subID']['sim']['num_sim_per_bin']

    # Generate the list of bin edges using the provided binning parameters
    bin_list = generate_bin_list(select_bin_min, select_bin_max,
                                 select_bin_interval)
    bin_list[0] += 1e-6
    bin_list[-1] -= 1e-6

    # Initialize an empty DataFrame with the same columns as df_subID
    df_expanded = pd.DataFrame(columns=df_subID.columns)

    # Mapping of windowIDs to corresponding indices
    largest_windowID = np.max(
        df_subID["windowID"].str.extract(r'win-(\d+)')[0].astype(int))
    windowID_map = {
        f"win-{str(i+1).zfill(2)}": i
        for i in range(0, largest_windowID)
    }

    # Convert "windowID" column to numeric for easy filtering
    df_subID["window_num"] = df_subID["windowID"].map(windowID_map)

    # Prepare lists for sim_sesID, sim_runID, sim_run_percent_censored, and sim_run_bin
    rm_run_th_list = []
    combo_list = []
    sim_sesID_list = []
    sim_runID_list = []
    sim_bld_list = []
    sim_run_percent_censored_list = []
    sim_run_FD_metric_list = []
    sim_run_DV_metric_list = []
    sim_run_bin_list = []

    # Loop through each % censored frames threshold
    for rm_run_th_value in sorted(combo_dict, reverse=True):
        # Loop through each unique combination
        for combo in combo_dict[rm_run_th_value]:
            # Define combo string value to be inserted
            # into output csv later
            combo_str = "l".join(combo)

            # Loop through each instance,
            # 1 instance = 1 simulation sesID
            for sim_ses_idx in range(num_sim_per_bin):
                # Define ses_key
                ses_key = f"ses-{str((sim_ses_idx + 1)).zfill(2)}"

                # Loop through each bin
                for sim_run_idx, bin_key in enumerate(combo):
                    # Define run_key
                    run_key = f"run-{str((sim_run_idx + 1)).zfill(2)}"

                    # Define sim_bld
                    sim_bld = f"{str((sim_run_idx + 1)).zfill(3)}"

                    # Extract selected instances
                    select_instance = select_dict[run_key][bin_key]

                    # If sufficient bin,
                    if select_instance is not None:
                        # Extracted selected window rows from window list
                        row_indices_for_run = select_instance[:, sim_ses_idx]
                        selected_rows = df_subID.loc[row_indices_for_run].copy(
                        )

                        # Select rows where "window_num" matches `row_indices_for_run`
                        selected_rows = df_subID[df_subID["window_num"].isin(
                            row_indices_for_run)].copy()

                        # Ensure the order matches `row_indices_for_run`
                        selected_rows["window_num"] = pd.Categorical(
                            selected_rows["window_num"],
                            categories=row_indices_for_run,
                            ordered=True)
                        selected_rows = selected_rows.sort_values("window_num")

                        # Append the selected rows to df_expanded
                        df_expanded = pd.concat([df_expanded, selected_rows],
                                                ignore_index=True)

                        # Compute FD & DV metrics for simulation run
                        sim_percent_censored = np.mean(
                            selected_rows["window_percent_censored"])
                        sim_FD_metric = np.mean(
                            selected_rows[f"window_mean_{FD_metric}"])
                        sim_DV_metric = np.mean(
                            selected_rows[f"window_mean_{DV_metric}"])

                        # Use pd.cut to assign each sim_run_percent_censored to a bin interval
                        bin_category = pd.cut([sim_percent_censored],
                                              bins=bin_list,
                                              right=True)

                        # Get the bin interval as a string
                        sim_run_bin = f"{int(bin_category[0].left)}_{int(bin_category[0].right)}"

                        # Append desired variables to lists
                        rm_run_th_list.extend([rm_run_th_value] *
                                              len(row_indices_for_run))
                        combo_list.extend([combo_str] *
                                          len(row_indices_for_run))
                        sim_sesID_list.extend([ses_key] *
                                              len(row_indices_for_run))
                        sim_runID_list.extend([run_key] *
                                              len(row_indices_for_run))
                        sim_bld_list.extend([sim_bld] *
                                            len(row_indices_for_run))
                        sim_run_percent_censored_list.extend(
                            [sim_percent_censored] * len(row_indices_for_run))
                        sim_run_FD_metric_list.extend([sim_FD_metric] *
                                                      len(row_indices_for_run))
                        sim_run_DV_metric_list.extend([sim_DV_metric] *
                                                      len(row_indices_for_run))
                        sim_run_bin_list.extend([sim_run_bin] *
                                                len(row_indices_for_run))

    # After the loop, append desired variables as new columns
    df_expanded['rm_run_th'] = rm_run_th_list
    df_expanded['combo'] = combo_list
    df_expanded['sim_sesID'] = sim_sesID_list
    df_expanded['sim_runID'] = sim_runID_list
    df_expanded['sim_bld'] = sim_bld_list
    df_expanded['sim_run_bin'] = sim_run_bin_list
    df_expanded['sim_run_percent_censored'] = sim_run_percent_censored_list
    df_expanded[f'sim_run_{FD_metric}'] = sim_run_FD_metric_list
    df_expanded[f'sim_run_{DV_metric}'] = sim_run_DV_metric_list

    return df_expanded


def compute_weighted_probability_vector(bin_min, bin_max, bin_interval,
                                        percent_censored_arr):
    """Compute a weighted probability vector for a set of percent censored values within a bin.

    Sub-bins `[bin_min, bin_max]` into `bin_interval`-wide bins, then assigns
    each element of `percent_censored_arr` an inverse-frequency weight based
    on which sub-bin it falls in (rarer sub-bins get higher weight), so that
    subsequent random selection favors a good spread of values across the
    sub-bins rather than clustering.

    Args:
        bin_min (int or float): Minimum value of the bin range.
        bin_max (int or float): Maximum value of the bin range.
        bin_interval (int or float): The interval used to create sub-bins.
        percent_censored_arr (np.ndarray): Array of percent censored values
            for which the probability vector is computed.

    Returns:
        np.ndarray: Weighted probability vector corresponding to the input censored values.
    """

    # Compute bin edges and digitize
    bin_edges = np.arange(bin_min, bin_max + bin_interval, bin_interval)
    bin_indices = np.digitize(percent_censored_arr, bin_edges, right=True)

    # Validate indices
    assert 0 not in bin_indices, f"ERROR: `bin_indices` contains elements <= `bin_min`"
    assert len(
        bin_edges
    ) not in bin_indices, f"ERROR: `bin_indices` contains elements > `bin_max`"

    # Shift indices down by 1 to avoid extra bin
    bin_indices = bin_indices - 1

    # Count occurrences of each bin and compute bin weights
    bin_counts = np.bincount(bin_indices, minlength=len(bin_edges) - 1)
    # If `bin_counts` > 0, set weights as 1/bin_counts
    bin_weights = np.divide(1, bin_counts, where=bin_counts > 0)
    # Else, set weights as 0
    bin_weights[bin_counts == 0] = 0

    # Normalize bin weights
    bin_weights /= bin_weights.sum()

    # Map bin weights to the corresponding bin indices
    p_vector = bin_weights[bin_indices]
    # Normalize the final probability vector
    p_vector /= p_vector.sum()

    return p_vector


def filter_zero_count_bin(pl_config_file, select_bin_counts,
                          sim_run_percent_censored):
    """Filter out selection bins with zero counts or an insufficient spread of % censored frames.

    Bins with 0 counts are dropped outright. For censored bins (excluding
    the 0% censored bin) with non-zero counts, the bin is also dropped if
    the actual range of % censored frames values that fall within it covers
    less than `keep_select_bin_percent_th`% of the bin's interval width
    (i.e. the samples are too tightly clustered to represent the full bin).

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        select_bin_counts (pd.Series): Bin counts (indexed by bin interval),
            e.g. from `compute_bin_count`.
        sim_run_percent_censored (np.ndarray): Mean % censored frames for
            each simulation run.

    Returns:
        pd.Series: Filtered bins with sufficient counts and spread.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    select_bin_interval = pl_config['select_subID']['sim'][
        'select_bin_interval']
    keep_select_bin_percent_th = pl_config['select_subID']['sim'][
        'keep_select_bin_percent_th']

    # Drop bins with 0 counts
    select_bin_counts = select_bin_counts[select_bin_counts > 0]
    sufficient_select_bins = select_bin_counts.copy()

    for bin_idx, bin_count in select_bin_counts.items():
        # If censored bin (i.e., if NOT 0% censored frames bin)
        if bin_idx.right != 0:
            bin_min, bin_max = bin_idx.left, bin_idx.right
            bin_mask = (sim_run_percent_censored
                        > bin_min) & (sim_run_percent_censored <= bin_max)
            sim_run_bin_percent_censored = sim_run_percent_censored[bin_mask]
            sim_run_bin_percent_censored_range = np.max(
                sim_run_bin_percent_censored) - np.min(
                    sim_run_bin_percent_censored)

            # If bin range does not have >= `keep_select_bin_percent_th`% of bin interval, discard this bin
            if (sim_run_bin_percent_censored_range /
                    select_bin_interval) * 100 < keep_select_bin_percent_th:
                sufficient_select_bins.drop(bin_idx, inplace=True)

    return sufficient_select_bins


def get_combination(pl_config_file, censored_bin_list):
    """Determine combinations of selection bins across simulation runs, for each rm_run_th threshold.

    For each % censored frames threshold (to discard a run), generates all
    combinations-with-replacement of `censored_bin_list` of length
    `num_sim_run_per_sim_ses` (i.e. the "bottom right triangle" of the grid
    of bin x bin x ... combinations, retaining combinations but not
    permutations). Each combination is filtered to only keep bins whose
    lower edge is below the threshold, and duplicate (order-independent)
    filtered combinations are dropped.

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        censored_bin_list (list): List of censored selection bin strings
            (e.g. "0_25"), excluding the non-censored bin.

    Returns:
        dict: Dictionary of combinations of bins and simulation runs, keyed
            by rm_run_th threshold value, each mapping to a list of
            bin-combination tuples.
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rm_run_th = pl_config['select_subID']['sim']['rm_run_th']
    num_sim_run_per_sim_ses = pl_config['select_subID']['sim'][
        'num_sim_run_per_sim_ses']

    # Init dictionary storing all final analysis
    combo_dict = {}
    largest_num_sim_run_per_sim_ses = 0

    for rm_run_th_idx, rm_run_th_value in enumerate(rm_run_th):
        # Init rm_run_th_value-level variables
        combo_dict[rm_run_th_value] = []
        seen_sets = set()

        # Obtain "bottom right triangle" (inclusive of diagonal)
        combinations = itertools.combinations_with_replacement(
            censored_bin_list, num_sim_run_per_sim_ses)

        # Loop through each combination to drop 1) combinations
        # that do not pass % censored frames threshold (to discard
        # run) & 2) repeated combinations.
        for combo in combinations:
            # Only keep simulation runs that do not exceed
            # % censored frames threshold (to discard run)
            pass_th_combo = [
                censored_bin for censored_bin in combo
                if int(censored_bin.split("_")[0]) < rm_run_th_value
            ]
            pass_th_combo_tuple = tuple(sorted(pass_th_combo))

            # If no exceed threshold & not seen before, append
            if pass_th_combo and (pass_th_combo_tuple not in seen_sets):
                combo_dict[rm_run_th_value].append(pass_th_combo)
                if len(pass_th_combo) > largest_num_sim_run_per_sim_ses:
                    largest_num_sim_run_per_sim_ses = len(pass_th_combo)

            # Append to seen set
            if pass_th_combo_tuple not in seen_sets:
                seen_sets.add(pass_th_combo_tuple)

    return combo_dict


def initialize_select_dict(pl_config_file, combo_dict):
    """Initialize an (empty) selection dictionary based on the provided bin combinations.

    Using the combinations for the loosest (max) rm_run_th threshold in
    `combo_dict`, builds a nested dictionary keyed by run (e.g. "run-01")
    then by bin (e.g. "0_25"), with each bin initialized to None to be
    filled in later by `select_censored_bin`/`select_noncensored_bin`.

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        combo_dict (dict): Dictionary (keyed by rm_run_th value) of lists of
            bin-combination tuples, as produced by `get_combination`.

    Returns:
        dict: Initialized dictionary with simulation run IDs and bin values,
            of the form {run_key: {bin_key: None}}.
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rm_run_th = pl_config['select_subID']['sim']['rm_run_th']

    # Initialize select_dict
    select_dict = {}

    # Get the loosest run threshold value (max value)
    loosest_rm_run_th_value = np.max(rm_run_th)

    # Loop through the combinations for the loosest run threshold
    for combo in combo_dict[loosest_rm_run_th_value]:
        for bin_idx, bin_value in enumerate(combo):
            # Create run key (e.g., "run-01", "run-02", etc.)
            run_key = f"run-{str((bin_idx + 1)).zfill(2)}"

            if run_key not in select_dict:
                select_dict[run_key] = {}

            # Create simulation run ID key for each bin value
            if bin_value not in select_dict[run_key]:
                select_dict[run_key][bin_value] = None

            # Sort simulation run IDs in order
            select_dict[run_key] = dict(sorted(select_dict[run_key].items()))

    return select_dict


def drop_non_unique_combination(combo_dict):
    """Drop non-unique bin combinations to avoid repeated analyses.

    Combinations are considered duplicates across thresholds too: iterates
    thresholds from loosest to strictest, and once a combination has been
    seen (at a looser threshold), it is dropped from all stricter
    thresholds.

    Args:
        combo_dict (dict): Dictionary (keyed by rm_run_th value) of lists of
            bin-combination tuples/lists, as produced by `get_combination`.

    Returns:
        dict: Updated dictionary with only unique bin combinations for each threshold.
    """
    seen_combinations = set()

    # Loop through discard run thresholds, from loosest to strictest thresholds.
    for key in sorted(combo_dict.keys(), reverse=True):
        unique_combinations = []

        # Loop through combinations
        for combination in combo_dict[key]:
            combo_tuple = tuple(combination)
            # If unique combination, append
            if combo_tuple not in seen_combinations:
                unique_combinations.append(combination)
                # Mark this combination as seen
                seen_combinations.add(combo_tuple)

        # Update `combo_dict`
        combo_dict[key] = unique_combinations

    return combo_dict


#############################################################################
# ===========================================================================
# SUB_ID LEVEL PLOT
# ===========================================================================
#############################################################################
def plot_1ses2run_no_rm_run_th_heatmap(pl_config_file, select_dict, dataset,
                                       subID):
    """Plot and save a heatmap of successfully-selected simulation run counts, for a single subID.

    Only applicable to the 2-simulation-runs-per-simulation-session set-up
    (i.e. no `rm_run_th` filtering, hence "no_rm_run_th" in the name). For
    every pair of motion bins (bin_1 for simulation run 1, bin_2 for
    simulation run 2), counts how many of the 2 runs had a successfully
    selected (non-None) instance for that bin in `select_dict`, and plots
    this count (0, 1, or 2) as a heatmap restricted to the upper-right
    triangle (bin_1 <= bin_2).

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        select_dict (dict): Nested dictionary of {run_key: {bin_key:
            selected window-index array or None}}, as filled in by
            `select_censored_bin`.
        dataset (str): Dataset name (used in the plot title/filename).
        subID (str): Subject ID (used in the plot title/filename).
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    select_bin_min = pl_config['select_subID']['sim']['select_bin_min']
    select_bin_max = pl_config['select_subID']['sim']['select_bin_max']
    select_bin_interval = pl_config['select_subID']['sim'][
        'select_bin_interval']
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']

    # Define output dir
    selected_heatmap_dir = os.path.join(rep_sel_subID_sim_dir, 'result',
                                        'selected', 'heatmap')
    os.makedirs(selected_heatmap_dir, exist_ok=True)

    # Define motion_bins
    motion_bins = [
        f"{i}_{i + select_bin_interval}"
        for i in range(select_bin_min, select_bin_max, select_bin_interval)
        if f"{i}_{i + select_bin_interval}"[0] != "-"
    ]
    num_bins = len(motion_bins)

    # Init heatmap array
    data = np.zeros((num_bins, num_bins))

    # Iterate through all combinations of simulation runs and motion bins
    for sim_run_1 in range(num_bins):
        for sim_run_2 in range(
                sim_run_1, num_bins
        ):  # Ensure we avoid repeated combinations in the top-left triangle
            bin_1 = motion_bins[sim_run_1]
            bin_2 = motion_bins[sim_run_2]

            # Count the non-None numpy arrays in select_dict for this bin combination
            count = 0

            # Check the presence of np arrays for each simulation run
            if select_dict["run-01"].get(bin_1) is not None:
                count += 1
            if select_dict["run-02"].get(bin_2) is not None:
                count += 1

            # Store the count in the data array
            data[sim_run_1, sim_run_2] = count

    # Now, we will plot the heatmap with the updated 'data' array
    # Create a mask for the top-left triangle (excluding diagonal)
    mask = np.tril(np.ones_like(data, dtype=bool), k=-1)

    # Create the heatmap
    plt.figure(figsize=(8, 6))
    heatmap = sns.heatmap(data,
                          cmap="Blues",
                          cbar=True,
                          xticklabels=motion_bins,
                          yticklabels=motion_bins,
                          linewidths=1,
                          linecolor='black',
                          square=True,
                          mask=mask,
                          vmin=0,
                          vmax=2,
                          cbar_kws={"ticks": [0, 1, 2]})

    # Modify the mask for top-left cells to be black (overwrites the default white)
    for i in range(num_bins):
        for j in range(i):
            plt.gca().add_patch(plt.Rectangle(
                (j, i), 1, 1,
                color='black'))  # Fill top-left triangle with black

    # Set title and labels
    plt.title(f"{dataset} {subID}", fontsize=16)
    plt.xlabel("Simulation Run 2", fontsize=12)
    plt.ylabel("Simulation Run 1", fontsize=12)

    # Invert the y-axis to have the first bin at the bottom
    plt.gca().invert_yaxis()

    # Get the color bar axis and label it
    cbar = heatmap.collections[0].colorbar
    cbar.set_label("Number of successfully simulated runs", fontsize=12)

    # Save heatmap
    heatmap_filename = os.path.join(selected_heatmap_dir,
                                    f"{dataset}_{subID}.png")
    plt.savefig(heatmap_filename, bbox_inches="tight")
    # plt.tight_layout()
    # plt.show()


#############################################################################
# ===========================================================================
# JOB SUBMISSION
# ===========================================================================
#############################################################################


def submit_select_subID_sim_job(timestamp, jobprefix, main_config_file,
                                pl_config_file, rep_sel_subID_sim_dir, dataset,
                                subID, std_dir, username):
    """Submit 1 job to run `select_subID_sim` for 1 subID.

    Deletes any pre-existing output files for this subID (candidate plot,
    selected csv/plot, rejected csv) before submitting, so that stale
    outputs from a previous run cannot be mistaken for fresh results.

    Args:
        timestamp (str): Timestamp string, currently unused within this
            function body (kept for a consistent job-submission function
            signature).
        jobprefix (str): Prefix used to build the job name.
        main_config_file (str): Path to the main config JSON file.
        pl_config_file (str): Path to the pipeline config JSON file.
        rep_sel_subID_sim_dir (str): Base directory for select_subID sim
            outputs and old-output cleanup.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        std_dir (str): Directory to save the job's stdout/stderr files to.
        username (str): Username, currently unused within this function body
            (kept for a consistent job-submission function signature).
    """

    # Job directories and paths
    os.makedirs(std_dir, exist_ok=True)

    # Delete old output(s)
    file_paths = [
        os.path.join(rep_sel_subID_sim_dir, 'result', 'candidate',
                     f"{dataset}_{subID}.png"),
        os.path.join(rep_sel_subID_sim_dir, 'result', 'selected',
                     f"{dataset}_{subID}.csv"),
        os.path.join(rep_sel_subID_sim_dir, 'result', 'selected',
                     f"{dataset}_{subID}.png"),
        os.path.join(rep_sel_subID_sim_dir, 'result', 'rejected',
                     f"{dataset}_{subID}.csv")
    ]
    for file_path in file_paths:
        if os.path.exists(file_path):
            print(f"Deleting existing file: {file_path}")
            os.remove(file_path)

    # Generate job-specific variables
    jobname = f"{jobprefix}_{dataset}_{subID}"
    stderr_file = os.path.join(std_dir, f"{jobname}.err")
    stdout_file = os.path.join(std_dir, f"{jobname}.out")

    # Command to call the filtering function
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; sys.path.append(\\\"{PROJECT_DIR}/select_subID/sim\\\"); "
        f"import CBIG_TMSmotion_select_subID_sim as sel_subID_sim; "
        f"sel_subID_sim.select_subID_sim(\\\"{main_config_file}\\\", \\\"{pl_config_file}\\\", "
        f"\\\"{dataset}\\\", \\\"{subID}\\\")'; "
        f"conda deactivate")

    # Submit job
    util.submit_job(walltime="24:00:00",
                    memory="50G",
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def submit_all_select_subID_sim_job(timestamp, jobprefix, main_config_file,
                                    pl_config_file, rep_sel_subID_sim_dir,
                                    list_file, std_dir, username):
    """Submit `select_subID_sim` jobs for every (dataset, subID) pair in a list file.

    Submits all jobs for a list of subIDs.
    Each job is for each subID.

    Args:
        timestamp (str): Timestamp string, passed through to
            `submit_select_subID_sim_job`.
        jobprefix (str): Prefix used to build each job's name.
        main_config_file (str): Path to the main config JSON file.
        pl_config_file (str): Path to the pipeline config JSON file.
        rep_sel_subID_sim_dir (str): Base directory for select_subID sim
            outputs and old-output cleanup.
        list_file (str): Path to a CSV file with 'dataset' and 'subID'
            columns (and a 'bld' column, read as str).
        std_dir (str): Directory to save each job's stdout/stderr files to.
        username (str): Username, passed through to
            `submit_select_subID_sim_job`.
    """

    # Submit job
    df = pd.read_csv(list_file, dtype={'bld': str})
    for (dataset, subID), _ in df.groupby(['dataset', 'subID']):
        submit_select_subID_sim_job(timestamp, jobprefix, main_config_file,
                                    pl_config_file, rep_sel_subID_sim_dir,
                                    dataset, subID, std_dir, username)


def check_select_subID_sim_job(timestamp, try_count, dataset, subID,
                               rep_sel_subID_sim_dir, all_job_success):
    """Check whether a single subID's `select_subID_sim` job completed successfully.

    A job is considered successful if its output csv exists in exactly one
    of the 'selected' or 'rejected' result directories (existing in both, or
    neither, is treated as a failure). On failure, appends an entry to a
    per-try failure-message file and a per-try failure-list file.

    Args:
        timestamp (str): Timestamp string identifying the output
            subdirectory to write failure files to.
        try_count (int): Current retry attempt number, used to name the
            failure-message/failure-list files.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        rep_sel_subID_sim_dir (str): Base directory for select_subID sim
            outputs.
        all_job_success (bool): Running success flag across all checked
            jobs so far; set to False if this subID's job failed.

    Returns:
        tuple: A 2-tuple of (all_job_success (bool): updated running success
            flag; fail_list (str): path to the per-try failure-list file).
    """
    # Define filename to check
    filename2check = f"{dataset}_{subID}.csv"

    # Paths to the directories
    timestamp_dir = os.path.join(rep_sel_subID_sim_dir, timestamp)
    os.makedirs(timestamp_dir, exist_ok=True)
    selected_dir = os.path.join(rep_sel_subID_sim_dir, 'result', 'selected',
                                'df')
    rejected_dir = os.path.join(rep_sel_subID_sim_dir, 'result', 'rejected',
                                'df')

    # Flags to track file presence
    in_selected = os.path.exists(os.path.join(selected_dir, filename2check))
    in_rejected = os.path.exists(os.path.join(rejected_dir, filename2check))

    # Check job output files for success
    id_success = True
    # Check if the file exists in exactly one directory
    if in_selected and in_rejected:
        # The file should not exist in both directories
        id_success = False
        all_job_success = False
        fail_msg = os.path.join(timestamp_dir, f"fail_msg_try_{try_count}.txt")
        # Handle failure message list
        if not os.path.exists(fail_msg):
            with open(fail_msg, 'w') as f:
                f.write("dataset,subID    error_message\n")
        with open(fail_msg, 'a') as f:
            f.write(
                f"{dataset},{subID},Error: {dataset}_{subID}.csv found in both 'selected' and 'rejected' directories\n"
            )
    elif not (in_selected or in_rejected):
        # The file must exist in at least one directory
        id_success = False
        all_job_success = False
        fail_msg = os.path.join(timestamp_dir, f"fail_msg_try_{try_count}.txt")
        # Handle failure message list
        if not os.path.exists(fail_msg):
            with open(fail_msg, 'w') as f:
                f.write("dataset,subID    error_message\n")
        with open(fail_msg, 'a') as f:
            f.write(
                f"{dataset},{subID},Error: {dataset}_{subID}.csv not found in 'selected' or 'rejected' directories\n"
            )

    # Return success status and path to the failure message file (if any)
    fail_list = os.path.join(timestamp_dir, f"fail_list_try_{try_count}.txt")
    if not id_success:
        # Handle failed ID list
        if not os.path.exists(fail_list):
            with open(fail_list, 'w') as f:
                f.write("dataset,subID\n")
        with open(fail_list, 'a') as f:
            f.write(f"{dataset},{subID}\n")

    return all_job_success, fail_list


def check_all_select_subID_sim_job(timestamp, try_count, list_file,
                                   rep_sel_subID_sim_dir):
    """Check `select_subID_sim` jobs for every (dataset, subID) pair in a list file.

    Args:
        timestamp (str): Timestamp string identifying the output
            subdirectory to write failure files to.
        try_count (int): Current retry attempt number, used to name the
            failure-message/failure-list files.
        list_file (str): Path to a CSV file with 'dataset' and 'subID'
            columns (and a 'bld' column, read as str).
        rep_sel_subID_sim_dir (str): Base directory for select_subID sim
            outputs.

    Returns:
        tuple: A 2-tuple of (all_job_success (bool): True if every job
            succeeded; fail_list (str or None): path to the per-try
            failure-list file, or None if no jobs were checked).
    """
    # Initialize success flag and failure list
    all_job_success = True
    fail_list = None

    # Iterate over each row in the DataFrame
    df = pd.read_csv(list_file, dtype={'bld': str})
    for (dataset, subID), _ in df.groupby(['dataset', 'subID']):
        # Call the check function for each dataset and subID
        all_job_success, fail_list = check_select_subID_sim_job(
            timestamp, try_count, dataset, subID, rep_sel_subID_sim_dir,
            all_job_success)

    return all_job_success, fail_list


#############################################################################
# ===========================================================================
# MULTIPLE SIMULATION RUNS IN SIMULATION SESSION
# ===========================================================================
#############################################################################


#############################################################################
# ===========================================================================
# GENERATE LISTS (SUMMARIZE POST-SIMULATION)
# ===========================================================================
#############################################################################
def generate_selected_rejected_subID_list_1(rep_sel_subID_sim_dir,
                                            out_dir,
                                            in_subID_list,
                                            choice='selected'):
    """Generate a list of subIDs that were selected or rejected during simulation.

    Generates a list for either selected or rejected subID.
    Each list comprises 2 columns: 'dataset' & 'subID'.
    Each row = each subID.

    Args:
        rep_sel_subID_sim_dir (str): Base directory for select_subID sim
            outputs, used to check for each subID's per-dataset result csv.
        out_dir (str): Output directory to save the generated list file to.
        in_subID_list (str): Path to a CSV file with 'dataset' and 'subID'
            columns (and 'bld'/'sim_bld' columns, read as str) listing
            candidate subIDs.
        choice (str, optional): Either 'selected' or 'rejected'; determines
            which result subdirectory is checked and the output filename.
            Defaults to 'selected'.

    Returns:
        str: Path to the generated '{choice}_subID_list_1.txt' file.
    """
    # Define the output file path
    os.makedirs(out_dir, exist_ok=True)
    output_file = os.path.join(out_dir, f'{choice}_subID_list_1.txt')

    # Overwrite any existing file by creating a new file with just the header
    with open(output_file, 'w') as f:
        f.write("dataset,subID\n")

    # Read the input CSV file
    in_df = pd.read_csv(in_subID_list, dtype={'bld': str, 'sim_bld': str})

    # Iterate through the grouped dataset and subID pairs
    for (dataset, subID), _ in in_df.groupby(['dataset', 'subID']):
        # Check if the path exists
        subID_path = os.path.join(rep_sel_subID_sim_dir, 'result', choice,
                                  'df', f"{dataset}_{subID}.csv")
        if os.path.exists(subID_path):
            # Append the dataset and subID to the output file
            with open(output_file, 'a') as f:
                f.write(f"{dataset},{subID}\n")
    return output_file


def stack_all_selected_subID_df(rep_sel_subID_sim_dir, in_subID_list):
    """Vertically stack the per-subID selected result dataframes into one dataframe.

    Args:
        rep_sel_subID_sim_dir (str): Base directory for select_subID sim
            outputs, used to look up each subID's selected result csv.
        in_subID_list (str): Path to a CSV file with 'dataset' and 'subID'
            columns (and 'rm_run_th'/'bld'/'sim_bld' columns, read as str)
            listing candidate subIDs.

    Returns:
        pd.DataFrame: Vertically stacked dataframe of all subIDs' selected
            results, or an empty DataFrame if none of the result files
            exist.
    """
    # Read the input CSV file
    in_df = pd.read_csv(in_subID_list,
                        dtype={
                            'rm_run_th': str,
                            'bld': str,
                            'sim_bld': str
                        })

    # Initialize an empty list to hold DataFrames
    df_list = []

    # Iterate through the grouped dataset and subID pairs
    for (dataset, subID), _ in in_df.groupby(['dataset', 'subID']):
        # Construct the file path for the subID CSV
        subID_path = os.path.join(rep_sel_subID_sim_dir, 'result', 'selected',
                                  'df', f"{dataset}_{subID}.csv")
        if os.path.exists(subID_path):
            # Read the CSV file and append to the list
            df_list.append(
                pd.read_csv(subID_path,
                            dtype={
                                'rm_run_th': str,
                                'bld': str,
                                'sim_bld': str
                            }))

    # Concatenate all DataFrames in the list
    if df_list:
        stacked_df = pd.concat(df_list, ignore_index=True)
        return stacked_df
    else:
        return pd.DataFrame()


def filter_out_combo_with_insufficient_bin(selected_subID_df):
    """Filter out (dataset, subID, rm_run_th, combo) groups with an insufficient number of runs.

    Filters out rows where the number of unique 'sim_bld' values does not match
    the expected number of runs based on the 'combo' column.

    Args:
        selected_subID_df (pd.DataFrame): Selected-subID window-level
            dataframe; must contain 'dataset', 'subID', 'rm_run_th', 'combo',
            and 'sim_bld' columns, where 'combo' bin names are joined by
            "l".

    Returns:
        pd.DataFrame: A filtered DataFrame with only valid subID groups.
    """
    filtered_rows = []

    # Group by dataset, subID, rm_run_th, and combo
    for (dataset, subID, rm_run_th,
         combo), group_df in selected_subID_df.groupby(
             ['dataset', 'subID', 'rm_run_th', 'combo']):
        # Compute the expected number of runs from the combo column
        expected_runs = combo.count(
            "l") + 1  # Number of bins = number of "l" + 1

        # Get the actual number of unique runs based on sim_bld
        actual_runs = group_df['sim_bld'].nunique()

        # Keep only if the expected number of runs match the actual runs
        if actual_runs == expected_runs:
            filtered_rows.append(group_df)

    # Concatenate the valid groups back into a DataFrame
    if filtered_rows:
        filtered_selected_subID_df = pd.concat(filtered_rows,
                                               ignore_index=False)
        filtered_selected_subID_df.reset_index(drop=True, inplace=True)
        return filtered_selected_subID_df

    # If no valid rows remain, return an empty DataFrame with the same columns
    return selected_subID_df.iloc[0:0]


def append_highest_sim_bld_column(selected_subID_df):
    """Add a 'highest_sim_bld' column recording the highest sim_bld value within each group.

    Adds a new column 'highest_sim_bld' to the DataFrame, which contains
    the highest 'sim_bld' value within each group of 'dataset', 'subID',
    'rm_run_th', 'combo', and 'sim_sesID'.

    Args:
        selected_subID_df (pd.DataFrame): Selected-subID window-level
            dataframe; must contain 'dataset', 'subID', 'rm_run_th', 'combo',
            'sim_sesID', and 'sim_bld' columns.

    Returns:
        pd.DataFrame: DataFrame with the 'highest_sim_bld' column added.
    """

    # Create 'highest_sim_bld' column
    selected_subID_df["highest_sim_bld"] = selected_subID_df.groupby(
        ["dataset", "subID", "rm_run_th", "combo",
         "sim_sesID"])["sim_bld"].transform("max")

    return selected_subID_df


def generate_selected_rejected_subID_list_2(selected_subID_df,
                                            selected_subID_list_1,
                                            rejected_subID_list_1, out_dir):
    """Update the selected/rejected subID lists after filtering leaves some subIDs with no bins.

    Updates the selected and rejected subID lists after filtering out subIDs no bins
    left for analysis.

    Any (dataset, subID) pair present in `selected_subID_list_1` but no
    longer present in `selected_subID_df` (e.g. because
    `filter_out_combo_with_insufficient_bin` dropped all of its rows) is
    moved from the selected list to the rejected list.

    Args:
        selected_subID_df (pd.DataFrame): Filtered selected-subID
            window-level dataframe; must contain 'dataset' and 'subID'
            columns.
        selected_subID_list_1 (str): Path to the original selected subID
            list csv (with 'dataset', 'subID' columns).
        rejected_subID_list_1 (str): Path to the original rejected subID
            list csv (with 'dataset', 'subID' columns).
        out_dir (str): Output directory to save the updated lists to.

    Returns:
        tuple: A 2-tuple of (selected_subID_list_2 (str): path to the
            updated selected subID list csv; rejected_subID_list_2 (str):
            path to the updated rejected subID list csv).
    """

    # Read the original selected and rejected subID lists
    selected_subID_list_1_df = pd.read_csv(selected_subID_list_1,
                                           dtype={
                                               'dataset': str,
                                               'subID': str
                                           })
    rejected_subID_list_1_df = pd.read_csv(rejected_subID_list_1,
                                           dtype={
                                               'dataset': str,
                                               'subID': str
                                           })

    # Get the unique (dataset, subID) pairs that remain after filtering
    remaining_subIDs = selected_subID_df[['dataset',
                                          'subID']].drop_duplicates()

    # Find subIDs that were in selected_subID_list_1 but are no longer in selected_subID_df
    removed_subIDs = selected_subID_list_1_df.merge(remaining_subIDs,
                                                    on=['dataset', 'subID'],
                                                    how='left',
                                                    indicator=True)
    removed_subIDs = removed_subIDs[removed_subIDs['_merge'] ==
                                    'left_only'].drop(columns=['_merge'])

    # Update selected_subID_list_2 by removing dropped subIDs
    selected_subID_list_2_df = selected_subID_list_1_df.merge(
        removed_subIDs, on=['dataset', 'subID'], how='left', indicator=True)
    selected_subID_list_2_df = selected_subID_list_2_df[
        selected_subID_list_2_df['_merge'] == 'left_only'].drop(
            columns=['_merge'])

    # Update rejected_subID_list_2 by adding the removed subIDs
    rejected_subID_list_2_df = pd.concat(
        [rejected_subID_list_1_df, removed_subIDs], ignore_index=True)

    # Define the output file path
    os.makedirs(out_dir, exist_ok=True)
    selected_subID_list_2 = os.path.join(out_dir, 'selected_subID_list_2.txt')
    rejected_subID_list_2 = os.path.join(out_dir, 'rejected_subID_list_2.txt')

    selected_subID_list_2_df.to_csv(selected_subID_list_2, index=False)
    rejected_subID_list_2_df.to_csv(rejected_subID_list_2, index=False)

    return selected_subID_list_2, rejected_subID_list_2


def generate_window_list_1(selected_subID_df, out_dir):
    """Save 1 window-list csv file per unique simulation run, mirroring the directory structure.

    Groups `selected_subID_df` by (dataset, subID, rm_run_th, combo,
    sim_run_bin, sim_sesID, sim_bld) so that each group is exactly the set
    of windows belonging to one simulation run, and saves each group to its
    own 'window_list_1.txt' under a
    out_dir/window_list_1/{dataset}/{subID}/{sim_run_bin}/{sim_sesID}/{sim_bld}/
    directory.

    Args:
        selected_subID_df (pd.DataFrame): Selected-subID window-level
            dataframe; must contain 'dataset', 'subID', 'rm_run_th', 'combo',
            'sim_run_bin', 'sim_sesID', and 'sim_bld' columns.
        out_dir (str): Base output directory under which the
            'window_list_1' subdirectory tree is created.
    """
    # Loop within each unique `sim_bld`.
    for _, group_df in selected_subID_df.groupby([
            'dataset', 'subID', 'rm_run_th', 'combo', 'sim_run_bin',
            'sim_sesID', 'sim_bld'
    ]):
        dataset = group_df['dataset'].iloc[0]
        subID = group_df['subID'].iloc[0]
        sim_run_bin = group_df['sim_run_bin'].iloc[0]
        sim_sesID = group_df['sim_sesID'].iloc[0]
        sim_bld = group_df['sim_bld'].iloc[0]

        # Define the output path for the specific combination
        out_path = os.path.join(out_dir, "window_list_1", dataset, subID,
                                sim_run_bin, sim_sesID, sim_bld)
        os.makedirs(out_path, exist_ok=True)

        # Save the DataFrame for the current combination
        out_file = os.path.join(out_path, "window_list_1.txt")
        group_df.to_csv(out_file, index=False)


def generate_sim_run_list_1(selected_subID_df, out_dir):
    """Save a run-level (1 row per simulation run) csv, deduplicated across rm_run_th/combo.

    Args:
        selected_subID_df (pd.DataFrame): Selected-subID window-level
            dataframe; must contain 'dataset', 'subID', 'sim_run_bin',
            'sim_sesID', and 'sim_bld' columns.
        out_dir (str): Output directory to save 'sim_run_list_1.txt' to.
    """
    # Drop windows within each unique `sim_bld`.
    # Drop repeated `sim_run_bin` across rm_run_ths & combos.
    unique_df = selected_subID_df.drop_duplicates(
        subset=['dataset', 'subID', 'sim_run_bin', 'sim_sesID', 'sim_bld'])

    # Save the deduplicated DataFrame
    os.makedirs(out_dir, exist_ok=True)
    sim_run_list_file = os.path.join(out_dir, "sim_run_list_1.txt")
    unique_df.to_csv(sim_run_list_file, index=False)


def generate_sim_run_list_2(selected_subID_df, out_dir,
                            num_sim_run_per_sim_ses):
    """Save a run-level csv like `generate_sim_run_list_1`, further filtered to needed sim_runIDs.

    Args:
        selected_subID_df (pd.DataFrame): Selected-subID window-level
            dataframe; must contain 'dataset', 'subID', 'sim_run_bin',
            'sim_sesID', and 'sim_bld' columns.
        out_dir (str): Output directory to save 'sim_run_list_2.txt' to.
        num_sim_run_per_sim_ses (int): Number of simulation runs per
            simulation session, passed to `util.filter_sim_runID` to drop
            currently-not-needed sim_runIDs.
    """
    # Drop windows within each unique `sim_bld`.
    # Drop repeated `sim_run_bin` across rm_run_ths & combos.
    unique_df = selected_subID_df.drop_duplicates(
        subset=['dataset', 'subID', 'sim_run_bin', 'sim_sesID', 'sim_bld'])

    # Drop currently-not-needed sim_runIDs
    unique_df = util.filter_sim_runID(unique_df, num_sim_run_per_sim_ses)

    # Save the deduplicated DataFrame
    os.makedirs(out_dir, exist_ok=True)
    sim_run_list_file = os.path.join(out_dir, "sim_run_list_2.txt")
    unique_df.to_csv(sim_run_list_file, index=False)


def generate_sim_run_list_3(selected_subID_df,
                            out_dir,
                            num_sim_run_per_sim_ses,
                            FD_metric,
                            DV_metric,
                            FD_metric_value,
                            DV_metric_value,
                            mode,
                            lite_flag=False):
    """Generate a run-level list expanded across varied FD/DV threshold combinations, and save to a file.

    Deduplicates to 1 row per simulation run (per rm_run_th/combo), filters
    to currently-needed sim_runIDs, then expands each remaining row into one
    row per (vary_FD_th, vary_DV_th) combination from
    `util.generate_FD_DV_combination`.

    Args:
        selected_subID_df (pd.DataFrame): Selected-subID window-level
            dataframe; must contain 'dataset', 'subID', 'rm_run_th', 'combo',
            'sim_run_bin', 'sim_sesID', and 'sim_bld' columns.
        out_dir (str): Output directory to save
            'sim_run_list_3_{mode}.txt' to.
        num_sim_run_per_sim_ses (int): Number of simulation runs per
            simulation session, passed to `util.filter_sim_runID`.
        FD_metric (str): Name of the FD metric (used to build the
            'vary_{FD_metric}_th' column name).
        DV_metric (str): Name of the DV metric (used to build the
            'vary_{DV_metric}_th' column name).
        FD_metric_value (list or array-like): Candidate FD threshold values,
            passed to `util.generate_FD_DV_combination`.
        DV_metric_value (list or array-like): Candidate DV threshold values,
            passed to `util.generate_FD_DV_combination`.
        mode (str): FD/DV combination mode (e.g. "slide" or "grid"), passed
            to `util.generate_FD_DV_combination` and used in the output
            filename.
        lite_flag (bool, optional): If True, restrict to only the least and
            most aggressive censoring thresholds (first and last of the
            full FD/DV combination grid) and sort rows for readability.
            Defaults to False.

    Returns:
        pd.DataFrame: The expanded, saved dataframe.
    """
    # Drop windows within each unique `sim_bld`
    # Keep repeated `sim_run_bin` across combos.
    unique_df = selected_subID_df.drop_duplicates(subset=[
        'dataset', 'subID', 'rm_run_th', 'combo', 'sim_run_bin', 'sim_sesID',
        'sim_bld'
    ])

    # Drop currently-not-needed sim_runIDs
    unique_df = util.filter_sim_runID(unique_df, num_sim_run_per_sim_ses)

    # Generate FD-DV combinations
    FD_DV_combination_list = util.generate_FD_DV_combination(
        FD_metric_value, DV_metric_value, mode)

    # For lite replication mode, only keep the least & most
    # aggressive censoring thresholds (first & last of the full grid)
    if lite_flag:
        FD_DV_combination_list = [
            FD_DV_combination_list[0], FD_DV_combination_list[-1]
        ]

    # Expand rows with FD and DV combinations
    expanded_rows = []
    for _, row in unique_df.iterrows():
        for vary_FD_th, vary_DV_th in FD_DV_combination_list:
            # Append FD and DV as new columns to each row
            row_copy = row.copy()
            row_copy[f'vary_{FD_metric}_th'] = vary_FD_th
            row_copy[f'vary_{DV_metric}_th'] = vary_DV_th
            expanded_rows.append(row_copy)

    # Create a new DataFrame with the expanded rows
    expanded_df = pd.DataFrame(expanded_rows)

    # Reorder columns
    first_column = [
        'dataset', 'subID', 'rm_run_th', 'combo', 'sim_run_bin', 'sim_sesID',
        'sim_runID', 'sim_bld', 'highest_sim_bld', f'vary_{FD_metric}_th',
        f'vary_{DV_metric}_th', f'sim_run_{FD_metric}', f'sim_run_{DV_metric}',
        'windowID', 'window_idx_start', 'window_idx_end', 'bld',
        'sim_run_percent_censored'
    ]
    expanded_df = util.reorder_df_columns(expanded_df, first_column)

    # For lite replication mode, sort rows for readability
    if lite_flag:
        sort_cols = [
            'dataset', 'subID', 'rm_run_th', 'combo', 'sim_run_bin',
            'sim_sesID', 'sim_runID', f'vary_{FD_metric}_th',
            f'vary_{DV_metric}_th'
        ]
        expanded_df = expanded_df.sort_values(by=sort_cols).reset_index(
            drop=True)

    # Save the expanded DataFrame
    os.makedirs(out_dir, exist_ok=True)
    sim_run_list_file = os.path.join(out_dir, f"sim_run_list_3_{mode}.txt")
    expanded_df.to_csv(sim_run_list_file, index=False)

    return expanded_df


def generate_selected_subID_non_groundtruth_df(selected_subID_df, out_dir):
    """Generate and save the list of non-ground-truth runs needed by the selected simulation windows.

    Processes the selected_subID_df to drop duplicate rows of repeated
    'dataset', 'subID', 'bld', & saves result to a file.

    Purpose of new list is to preprocess non-ground truth runs
    that are required by the selected simulation run's windows.
    (and not waste compute resources preprocessing non-needed
    non-ground truth runs)

    Args:
        selected_subID_df (pd.DataFrame): Selected-subID window-level
            dataframe; must contain 'dataset', 'subID', and 'bld' columns.
        out_dir (str): Output directory to save 'ngt_list_2.txt' to.

    Returns:
        pd.DataFrame: The deduplicated, saved dataframe.
    """
    # Drop duplicates based on specified columns
    ngt_2_df = selected_subID_df.drop_duplicates(
        subset=['dataset', 'subID', 'bld'])

    # Ensure the output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # Reorder columns
    first_column = ['dataset', 'subID', 'bld']
    ngt_2_df = util.reorder_df_columns(ngt_2_df, first_column)

    # Save the DataFrame to a file in the output directory
    output_path = os.path.join(out_dir, "ngt_list_2.txt")
    ngt_2_df.to_csv(output_path, index=False)

    print(f"ngt_list_2.txt saved to {output_path}")

    return ngt_2_df


def filter_out_sim_run_bin(df, select_bin_min, select_bin_interval):
    """Filter out rows belonging to the lowest (non-censored) selection bin.

    Args:
        df (pd.DataFrame): Dataframe containing a 'sim_run_bin' column
            (e.g. "0_25").
        select_bin_min (int): Minimum selection bin edge (e.g. 0), used to
            build the target bin string to filter out.
        select_bin_interval (int): Selection bin interval (width), used to
            build the target bin string to filter out.

    Returns:
        pd.DataFrame: `df` with rows belonging to the
            '{select_bin_min}_{select_bin_min + select_bin_interval}' bin
            removed.
    """
    # Create the target bin string
    target_bin = f"{select_bin_min}_{select_bin_min + select_bin_interval}"

    # Filter out rows where sim_run_bin matches the target bin
    filtered_df = df[df['sim_run_bin'] != target_bin]

    return filtered_df


#############################################################################
# ===========================================================================
# GENERATE PLOTS (SUMMARIZE POST-SIMULATION)
# ===========================================================================
#############################################################################
def plot_FD_DV_percent_censored_hexbin(selected_subID_df, out_dir, FD_metric,
                                       DV_metric):
    """Plot and save hexbin plots of FD/DV metric vs. % censored frames across simulation runs.

    Generates 2 hexbin plots (1 for FD_metric, 1 for DV_metric).
    Each hexbin plot will show the relationship between the FD/DV_metric
    and % censored frames (within sim_runID).

    We expect the higher the FD/DV_metric, the higher the % censored
    frames.

    Args:
        selected_subID_df (pd.DataFrame): Selected-subID window-level
            dataframe; must contain 'dataset', 'subID', 'sim_run_bin',
            'sim_sesID', 'sim_runID', 'sim_run_percent_censored', and
            'window_mean_{FD_metric}'/'window_mean_{DV_metric}' columns.
        out_dir (str): Output directory to save the hexbin plot images to.
        FD_metric (str): Name of the FD metric (used to build the
            'window_mean_{FD_metric}' column name and plot labels).
        DV_metric (str): Name of the DV metric (used to build the
            'window_mean_{DV_metric}' column name and plot labels).
    """
    # Drop duplicates to ensure unique simulation runs
    unique_runs_df = selected_subID_df.drop_duplicates(
        subset=['dataset', 'subID', 'sim_run_bin', 'sim_sesID', 'sim_runID'])

    # Sort and reset index for consistency
    unique_runs_df = unique_runs_df.sort_values(by=['dataset', 'subID'])
    unique_runs_df.reset_index(drop=True, inplace=True)

    # Metrics to evaluate
    metrics = [FD_metric, DV_metric]

    # Calculate the sim_run_mean_{metric} for each subset
    for metric in metrics:
        # Create the metric column (e.g., 'window_mean_FD', 'window_mean_DVARS', etc.)
        metric_column = f'window_mean_{metric}'

        # Group by the relevant columns and calculate the mean for each group
        means_df = unique_runs_df.groupby(
            ['dataset', 'subID', 'sim_run_bin', 'sim_sesID',
             'sim_runID'])[metric_column].mean().reset_index()

        # Merge back to get the mean values in the original dataframe
        unique_runs_df = pd.merge(
            unique_runs_df,
            means_df,
            on=['dataset', 'subID', 'sim_run_bin', 'sim_sesID', 'sim_runID'],
            suffixes=('', f'_mean_{metric}'))

        # Plot hexbin for each metric
        plt.figure(figsize=(8, 6))
        plt.hexbin(unique_runs_df['sim_run_percent_censored'],
                   unique_runs_df[f'{metric_column}_mean_{metric}'],
                   gridsize=15,
                   cmap='Greens')
        plt.colorbar(label='Counts')
        plt.xlabel(
            'Mean % Censored Frames (within each selected simulated run)')
        plt.ylabel(f'Mean {metric} (within each selected simulated run)')
        plt.title(f'Distribution of mean % censored frames & {metric}')

        # Save the plot
        os.makedirs(out_dir, exist_ok=True)
        output_path = os.path.join(out_dir,
                                   f'percent_censored_{metric}_dist.png')
        plt.savefig(output_path, dpi=300)
        # plt.show()
        # plt.close()


def plot_violin_percent_censored(selected_subID_df, out_dir, FD_metric,
                                 DV_metric, bin_max, bin_interval):
    """Plot and save a violin plot of mean % censored frames distribution, one violin per subID.

    Generates violin plot to describe mean % censored frames
    (within sim_runID) distribution for each selected subID.

    Args:
        selected_subID_df (pd.DataFrame): Selected-subID window-level
            dataframe; must contain 'dataset', 'subID', and
            'sim_run_percent_censored' columns.
        out_dir (str): Output directory to save the violin plot image to.
        FD_metric (str): Name of the FD metric; only used to decide the
            plot height (taller if not 'FD').
        DV_metric (str): Name of the DV metric, currently unused within this
            function body.
        bin_max (int): Maximum y-axis (% censored frames) tick value.
        bin_interval (int): Interval between y-axis ticks.
    """

    # Set Arial font
    plt.rcParams['font.family'] = 'Arial'

    # Ensure the dataset and subID columns are available for grouping
    if 'dataset' not in selected_subID_df.columns or 'subID' not in selected_subID_df.columns:
        raise ValueError(
            "The dataframe must contain 'dataset' and 'subID' columns.")

    # Create a column for the plot label (combining dataset and subID)
    selected_subID_df['plot_label'] = selected_subID_df[
        'dataset'] + " " + selected_subID_df['subID']

    # Create the violin plot
    height = 8 if FD_metric == 'FD' else 9
    plt.figure(figsize=(12, height))
    sns.violinplot(
        data=selected_subID_df,
        x='plot_label',
        y='sim_run_percent_censored',
        inner='point',  # Plot points for each unique simulated run
        scale='area',  # Scale the area of each violin by count
        bw='silverman',  # Bandwidth estimator for smoothing
        cut=0)  # No extrapolation beyond the extreme data points

    # Customize the plot's title and axis labels
    plt.xlabel('Participant', fontsize=25, labelpad=17.5)
    plt.ylabel('Mean % Censored Frames', fontsize=25, labelpad=10)

    # Customize y-ticks
    plt.yticks(range(0, bin_max + 10, bin_interval), fontsize=22)

    # Rotate x-axis labels to avoid overlap
    plt.xticks(rotation=65, ha='right', rotation_mode='anchor', fontsize=22)

    # Add gridlines to the y-axis for readability
    plt.grid(axis='y')

    # Adjust layout to prevent clipping of labels
    plt.tight_layout()

    # Save the plot to the specified output directory
    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, 'percent_censored_violin.png')
    plt.savefig(output_path, dpi=300)
    plt.close()
    # plt.show()


def plot_1ses2run_no_rm_run_th_group_heatmap(pl_config_file,
                                             selected_subID_list):
    """Plot and save a group-level heatmap of selected-participant counts across motion bin pairs.

    Generates 4x4 grid showing number of selected participants.

    Only applicable to the 2-simulation-runs-per-simulation-session set-up
    (i.e. no `rm_run_th` filtering, hence "no_rm_run_th" in the name). Loads
    each subID's `select_dict` pickle (saved by `select_subID_sim`), and for
    every pair of motion bins (bin_1 for simulation run 1, bin_2 for
    simulation run 2), counts across all subIDs how many had both bins
    successfully selected (non-None) in their `select_dict`, then plots this
    count as an annotated heatmap restricted to the upper-right triangle
    (bin_1 <= bin_2).

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
        selected_subID_list (str): Path to a CSV file with 'dataset' and
            'subID' columns (and 'bld'/'sim_bld' columns, read as str)
            listing the final selected subIDs.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    select_bin_min = pl_config['select_subID']['sim']['select_bin_min']
    select_bin_max = pl_config['select_subID']['sim']['select_bin_max']
    select_bin_interval = pl_config['select_subID']['sim'][
        'select_bin_interval']
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']

    # Set Arial font
    plt.rcParams['font.family'] = 'Arial'

    # Read the selected_subID_list CSV
    in_df = pd.read_csv(selected_subID_list,
                        dtype={
                            'bld': str,
                            'sim_bld': str
                        })

    # Collect valid .pkl file paths
    pkl_files = []
    for (dataset, subID), _ in in_df.groupby(['dataset', 'subID']):
        # Construct the file path for the subID .pkl file
        subID_path = os.path.join(rep_sel_subID_sim_dir, 'result', 'selected',
                                  'select_dict', f"{dataset}_{subID}.pkl")
        if os.path.exists(subID_path):
            pkl_files.append(subID_path)

    if not pkl_files:
        print("No valid .pkl files found for selected subIDs.")

    # Define output directory for group heatmap
    output_dir = os.path.join(rep_sel_subID_sim_dir, "result", "summary")
    os.makedirs(output_dir, exist_ok=True)

    # Define motion bins
    motion_bins = [
        f"{i}_{i + select_bin_interval}"
        for i in range(select_bin_min, select_bin_max, select_bin_interval)
        if f"{i}_{i + select_bin_interval}"[0] != "-"
    ]
    num_bins = len(motion_bins)

    # Init heatmap array for group-level data
    data = np.zeros((num_bins, num_bins))

    # Iterate through all subID `select_dict` files and aggregate data
    # Update the heatmap data
    for sim_run_1 in range(num_bins):
        for sim_run_2 in range(sim_run_1, num_bins):
            bin_1 = motion_bins[sim_run_1]
            bin_2 = motion_bins[sim_run_2]

            for pkl_file in pkl_files:
                with open(pkl_file, "rb") as f:
                    select_dict = pickle.load(f)

                    # Count the non-None numpy arrays in select_dict for this bin combination
                    if select_dict["run-01"].get(
                            bin_1) is not None and select_dict["run-02"].get(
                                bin_2) is not None:
                        data[sim_run_1, sim_run_2] += 1
    data = data.T

    # Create a mask for the top-left triangle (excluding diagonal)
    mask = np.triu(np.ones_like(data, dtype=bool), k=1)

    formatted_lmotion_bins = [
        f"{bin.replace('_', '-') + '%'}" if 'l' not in bin else
        f"{bin.split('l')[0].replace('_', '-') + '%'} | {bin.split('l')[1].replace('_', '-') + '%'}"
        for bin in motion_bins
    ]

    # Create the heatmap
    plt.figure(figsize=(8, 6))
    heatmap = sns.heatmap(data,
                          cmap="Blues",
                          cbar=True,
                          xticklabels=formatted_lmotion_bins,
                          yticklabels=formatted_lmotion_bins,
                          linewidths=1,
                          linecolor='black',
                          square=True,
                          mask=mask,
                          vmin=0,
                          vmax=data.max(),
                          cbar_kws={
                              "orientation": "vertical",
                              "pad": 0.05,
                              "shrink": 0.9
                          },
                          annot=True,
                          annot_kws={"size": 14})

    # Modify the mask for top-left cells to be black (overwrites the default white)
    for i in range(num_bins):
        for j in range(i):
            plt.gca().add_patch(plt.Rectangle(
                (i, j), 1, 1,
                color='black'))  # Fill top-left triangle with black

    # Set title and labels
    plt.xlabel("Simulation Run 2: % Censored Frames", fontsize=17, labelpad=10)
    plt.ylabel("Simulation Run 1: % Censored Frames", fontsize=17, labelpad=10)

    # Modify colorbar
    cbar = heatmap.collections[0].colorbar
    mappable = heatmap.collections[0]
    cbar.set_ticks([])
    vmin, vmax = mappable.get_clim()
    cbar_text_shift = 0.02  # shift slightly outside
    cbar.ax.text(0.5,
                 0 - cbar_text_shift,
                 f"{vmin:.0f}",
                 ha='center',
                 va='top',
                 transform=cbar.ax.transAxes,
                 fontsize=15)
    cbar.ax.text(0.5,
                 1 + cbar_text_shift,
                 f"{vmax:.0f}",
                 ha='center',
                 va='bottom',
                 transform=cbar.ax.transAxes,
                 fontsize=15)
    cbar.outline.set_edgecolor("black")
    cbar.outline.set_linewidth(1.5)
    cbar.set_label("Number of Selected Participants", fontsize=15, labelpad=10)

    # 1. Move x-axis labels and title to top
    heatmap.xaxis.set_ticks_position('top')  # Move ticks to top
    heatmap.xaxis.set_label_position('top')  # Move label to top

    # 2. Remove tick marks sticking out
    heatmap.tick_params(axis='both', which='both', length=0, labelsize=16)

    # Save heatmap
    heatmap_filename = os.path.join(output_dir,
                                    "selected_subID_group_heatmap.png")
    plt.savefig(heatmap_filename, bbox_inches="tight")
    plt.close()
    # plt.show()


#############################################################################
# ===========================================================================
# GENERATE TXT SUMMARY REPORT (SUMMARIZE POST-SIMULATION)
# ===========================================================================
#############################################################################
def write_summary_report(pass_gt_subID_list, selected_subID_list,
                         rejected_subID_list, FD_metric, DV_metric, out_dir):
    """Write a text summary report of subID counts at each ground truth / simulation stage.

    Generates summary report of number of subIDs that passed/failed
    different ground truth/simulation thresholds.

    Args:
        pass_gt_subID_list (str): Path to a CSV file listing subIDs that
            passed the ground truth threshold (with 'bld'/'sim_bld' columns,
            read as str).
        selected_subID_list (str): Path to a CSV file listing the final
            selected subIDs.
        rejected_subID_list (str): Path to a CSV file listing subIDs
            rejected by the simulation threshold.
        FD_metric (str): Name of the FD metric, currently unused within this
            function body.
        DV_metric (str): Name of the DV metric, currently unused within this
            function body.
        out_dir (str): Output directory to save 'summary_report.txt' to.
    """
    # Read input files
    pass_gt_subID_df = pd.read_csv(pass_gt_subID_list,
                                   dtype={
                                       'bld': str,
                                       'sim_bld': str
                                   })
    selected_subID_df = pd.read_csv(selected_subID_list,
                                    dtype={
                                        'bld': str,
                                        'sim_bld': str
                                    })
    rejected_subID_df = pd.read_csv(rejected_subID_list,
                                    dtype={
                                        'bld': str,
                                        'sim_bld': str
                                    })

    # Prepare the output path for the summary report
    summary_report_path = os.path.join(out_dir, 'summary_report.txt')

    with open(summary_report_path, 'w') as file:
        # Record retention stats
        file.write('# ====================================\n')
        file.write('# SAMPLE COUNT:\n')
        file.write('# ====================================\n')
        file.write(
            f'No. subIDs passed ground truth threshold: {len(pass_gt_subID_df)}\n'
        )
        file.write(
            f'No. subIDs failed simulation threshold: {len(rejected_subID_df)}\n'
        )
        file.write(f'No. final selected subIDs: {len(selected_subID_df)}\n')
        file.write('\n')


def generate_select_subID_summary(pl_config_file):
    """Wrapper that post-processes simulation results and generates the summary report and plots.

    This is a wrapper function that generates summary report &
    plots for selected subIDs.

    Stacks all per-subID selected result dataframes, filters out
    insufficient-bin combinations, computes the `highest_sim_bld` column,
    optionally restricts to the 2-simulation-run "pairs"/"dropped-1-run"
    subsets and to a "lite" replication-mode subset, then re-derives the
    selected/rejected subID lists and generates the various downstream
    window/run lists, plots, and the text summary report.

    Args:
        pl_config_file (str): Path to the pipeline config JSON file.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rep_sel_subID_ngt_dir = pl_config['path']['rep_sel_subID_ngt_dir']
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    select_bin_min = pl_config['select_subID']['sim']['select_bin_min']
    select_bin_max = pl_config['select_subID']['sim']['select_bin_max']
    select_bin_interval = pl_config['select_subID']['sim'][
        'select_bin_interval']
    num_sim_run_per_sim_ses = pl_config['select_subID']['sim'][
        'num_sim_run_per_sim_ses']
    FD_metric_value = pl_config['TMS']['FD_metric_value']
    DV_metric_value = pl_config['TMS']['DV_metric_value']
    lite_flag = pl_config['pipeline'].get('mode', 'full') == 'lite'

    # Assign input files
    subID_list = os.path.join(rep_sel_subID_ngt_dir, 'subID_list_1.txt')

    # Define output dir
    out_dir = os.path.join(rep_sel_subID_sim_dir, 'result', 'summary')

    # Generate selected and rejected subject ID lists
    selected_subID_list_1 = generate_selected_rejected_subID_list_1(
        rep_sel_subID_sim_dir, out_dir, subID_list, choice='selected')
    rejected_subID_list_1 = generate_selected_rejected_subID_list_1(
        rep_sel_subID_sim_dir, out_dir, subID_list, choice='rejected')

    # Stack all selected subject ID dataframes
    selected_subID_df_1 = stack_all_selected_subID_df(rep_sel_subID_sim_dir,
                                                      subID_list)

    # Filter out combos that do not have expected number of bins
    selected_subID_df_2 = filter_out_combo_with_insufficient_bin(
        selected_subID_df_1)

    # Create `highest_sim_bld` column
    selected_subID_df_2 = append_highest_sim_bld_column(selected_subID_df_2)

    # Keep subID's bins for both 1) bins with only both simulation run pairs,
    # and 2)"dropped-1-simulation-run" analysis.
    # For "dropped-1-simulation-run" analysis,
    # these subID's bins need to fulfill 2 conditions:
    # 1) SubID originally needed to have 60 pairs of simulation runs.
    # 2) These bins need to be from those that discard a simulation run
    #    due to exceeding the `rm_run_th` threshold of say, 50%.
    if num_sim_run_per_sim_ses == 2:
        # Keep bins with only both simulation run pairs
        selected_subID_df_2_pairs = selected_subID_df_2[
            selected_subID_df_2['highest_sim_bld'] == "002"].reset_index(
                drop=True)

        # Keep subID's bins for "dropped-1-simulation-run" analysis.
        # Step 1: Identify dataset and subID pairs that meet both conditions
        valid_combos = {
            "0_25l50_75", "25_50l50_75", "0_25l75_100", "25_50l75_100"
        }
        valid_pairs = selected_subID_df_2.loc[
            (selected_subID_df_2['highest_sim_bld'] == "002") &
            (selected_subID_df_2['combo'].isin(valid_combos)),
            ['dataset', 'subID']].drop_duplicates()
        # Step 2: Keep only rows where dataset and subID are in valid_pairs
        selected_subID_df_2_dropped = selected_subID_df_2.merge(
            valid_pairs, on=['dataset', 'subID'])
        # Step 3: Keep only bins with 1 remaining simulation run after filtering by `rm_run_th`
        selected_subID_df_2_dropped = selected_subID_df_2_dropped[
            selected_subID_df_2_dropped['rm_run_th'] == "50"].reset_index(
                drop=True)

        # Concatenate both dataframes
        selected_subID_df_2 = pd.concat(
            [selected_subID_df_2_pairs, selected_subID_df_2_dropped])

    # For lite replication mode: restrict to 1 designated motion bin `combo`
    # per subID, then the first `lite_num_sim_ses` simulated sessions. This
    # must stay a post-hoc restriction here (not baked into the earlier raw
    # window-selection step, e.g. via a smaller `num_sim_per_bin`), so that
    # everything kept is identical to what `full` mode would also produce
    # for those same combo(s)/session(s).
    if lite_flag:
        lite_combo = pl_config['select_subID']['sim']['lite_combo']
        subID_to_lite_combo = pd.Series(lite_combo)
        selected_subID_df_2 = selected_subID_df_2[
            selected_subID_df_2['combo'] == selected_subID_df_2['subID'].map(
                subID_to_lite_combo)].reset_index(drop=True)

        lite_num_sim_ses = pl_config['select_subID']['sim']['lite_num_sim_ses']
        keep_sim_sesID = sorted(
            selected_subID_df_2['sim_sesID'].unique())[:lite_num_sim_ses]
        selected_subID_df_2 = selected_subID_df_2[
            selected_subID_df_2['sim_sesID'].isin(keep_sim_sesID)].reset_index(
                drop=True)

    # Save `selected_subID_df_2` for debugging purposes
    selected_subID_df_2.to_csv(os.path.join(out_dir,
                                            'selected_subID_df_2.csv'),
                               index=False)

    # After `filter_out_combo_with_insufficient_bin`, some
    # subIDs might be dropped. Update selected_subID_list &
    # rejected_subID_list.
    selected_subID_list_2, rejected_subID_list_2 = generate_selected_rejected_subID_list_2(
        selected_subID_df_2, selected_subID_list_1, rejected_subID_list_1,
        out_dir)

    # Save 1 csv file per simulation run (each row = each window)
    generate_window_list_1(selected_subID_df_2, out_dir)

    # Save 1 csv file (each row = each simulation run)
    generate_sim_run_list_1(selected_subID_df_2, out_dir)

    # Save same as sim_run_list_1 except filtered out currently-not-needed sesIDs
    generate_sim_run_list_2(selected_subID_df_2, out_dir,
                            num_sim_run_per_sim_ses)

    # Save same as sim_run_list_1 except filtered out currently-not-needed sesIDs,
    # and contains multiple simulation run rows for each different varied FD metric
    # & DV metric thresholds
    generate_sim_run_list_3(selected_subID_df_2,
                            out_dir,
                            num_sim_run_per_sim_ses,
                            FD_metric,
                            DV_metric,
                            FD_metric_value,
                            DV_metric_value,
                            mode="slide",
                            lite_flag=lite_flag)
    generate_sim_run_list_3(selected_subID_df_2,
                            out_dir,
                            num_sim_run_per_sim_ses,
                            FD_metric,
                            DV_metric,
                            FD_metric_value,
                            DV_metric_value,
                            mode="grid",
                            lite_flag=lite_flag)

    # Generate non-groundtruth DataFrame and save it
    generate_selected_subID_non_groundtruth_df(selected_subID_df_2, out_dir)

    # Plot violin plots for percent censored
    plot_violin_percent_censored(selected_subID_df_2, out_dir, FD_metric,
                                 DV_metric, select_bin_max,
                                 select_bin_interval)

    # Filter out non-censored simulation run bins
    # for hexbin plot
    filtered_selected_subID_df = filter_out_sim_run_bin(
        selected_subID_df_2, select_bin_min, select_bin_interval)

    # Plot hexbin for FD and DV percent censored
    plot_FD_DV_percent_censored_hexbin(filtered_selected_subID_df, out_dir,
                                       FD_metric, DV_metric)

    # Plot group selected subID heatmap for 2 simulation runs
    # within 1 simulation session set-up
    if num_sim_run_per_sim_ses == 2:
        plot_1ses2run_no_rm_run_th_group_heatmap(pl_config_file,
                                                 selected_subID_list_2)

    # Write summary report
    write_summary_report(subID_list, selected_subID_list_2,
                         rejected_subID_list_2, FD_metric, DV_metric, out_dir)

    print("Select_subID report generation complete. Outputs saved to:",
          out_dir)


def submit_generate_select_subID_summary_job(timestamp, jobprefix,
                                             pl_config_file,
                                             rep_sel_subID_sim_dir, std_dir,
                                             username):
    """Submit 1 job to run `generate_select_subID_summary`.

    Deletes (restarts) the existing 'result/summary' output directory before
    submitting, so that stale outputs from a previous run cannot be mistaken
    for fresh results.

    Args:
        timestamp (str): Timestamp string, currently unused within this
            function body (kept for a consistent job-submission function
            signature).
        jobprefix (str): Used directly as the job name.
        pl_config_file (str): Path to the pipeline config JSON file.
        rep_sel_subID_sim_dir (str): Base directory for select_subID sim
            outputs; its 'result/summary' subdirectory is reset.
        std_dir (str): Directory to save the job's stdout/stderr files to.
        username (str): Username, currently unused within this function body
            (kept for a consistent job-submission function signature).
    """
    # Job directories and paths
    os.makedirs(std_dir, exist_ok=True)

    # Delete old output(s)
    summary_dir = os.path.join(rep_sel_subID_sim_dir, 'result', 'summary')
    util.restart_dir(summary_dir)

    # Generate job-specific variables
    jobname = f"{jobprefix}"
    stderr_file = os.path.join(std_dir, f"{jobname}.err")
    stdout_file = os.path.join(std_dir, f"{jobname}.out")

    # Command to call the filtering function
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; sys.path.append(\\\"{PROJECT_DIR}/select_subID/sim\\\"); "
        f"import CBIG_TMSmotion_select_subID_sim as sel_subID_sim; "
        f"sel_subID_sim.generate_select_subID_summary(\\\"{pl_config_file}\\\")'; "
        f"conda deactivate")

    # Submit job
    util.submit_job(walltime="1:00:00",
                    memory="8G",
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def check_generate_select_subID_summary_job(timestamp, try_count,
                                            rep_sel_subID_sim_dir):
    """Check whether the `generate_select_subID_summary` job completed successfully.

    Checks `generate_select_subID_summary` job output files to
    verify if job is completed successfully.

    Verifies that a fixed list of expected output filenames all exist under
    the 'result/summary' directory. On any missing file, appends an entry to
    a per-try failure-message file.

    Args:
        timestamp (str): Timestamp string identifying the output
            subdirectory to write the failure-message file to.
        try_count (int): Current retry attempt number, used to name the
            failure-message file.
        rep_sel_subID_sim_dir (str): Base directory for select_subID sim
            outputs.

    Returns:
        bool: True if all expected output files exist, False otherwise.
    """
    # Assume job successful
    all_job_success = True

    # Define file names to check if exist
    filename2check = [
        "ngt_list_2.txt", "percent_censored_violin.png",
        "selected_subID_list_1.txt", "rejected_subID_list_1.txt",
        "summary_report.txt", "sim_run_list_1.txt", "sim_run_list_2.txt",
        "sim_run_list_3_slide.txt"
    ]

    # Define output dir
    timestamp_dir = os.path.join(rep_sel_subID_sim_dir, timestamp)
    os.makedirs(timestamp_dir, exist_ok=True)

    # Define output file
    fail_msg = os.path.join(timestamp_dir, f"fail_msg_try_{try_count}.txt")

    # Check job output files for success
    for filename in filename2check:
        exist = util.check_job_output(
            f"{os.path.join(rep_sel_subID_sim_dir, 'result', 'summary', filename)}"
        )
        if not exist:
            all_job_success = False
            # Handle failure message list
            if not os.path.exists(fail_msg):
                with open(fail_msg, 'w') as f:
                    f.write("error_message\n")
            with open(fail_msg, 'a') as f:
                f.write(f"{filename}\n")

    return all_job_success
