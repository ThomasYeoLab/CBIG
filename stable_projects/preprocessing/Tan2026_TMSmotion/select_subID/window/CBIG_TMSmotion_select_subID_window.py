'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
"""
Examples:

python $TAN2026_TMSMOTION_DIR/select_subID/window/CBIG_TMSmotion_select_subID_window.py \
    $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json

python $TAN2026_TMSMOTION_DIR/select_subID/window/CBIG_TMSmotion_select_subID_window.py \
    $TAN2026_TMSMOTION_DIR/replication/pl2/config/pipeline.json

"""
# ===========================================================================
# IMPORTS
# ===========================================================================
import os
import sys
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

TAN2026_TMSMOTION_DIR = os.environ.get('TAN2026_TMSMOTION_DIR')
sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'util'))
import CBIG_TMSmotion_util as util


# ===========================================================================
# FUNCTIONS
# ===========================================================================
def generate_window_list(main_config_file, pl_config_file, preproc_type,
                         run_type):
    """Generate a list of sliding-window motion QC windows for each subID.

    Reads the ground-truth-passing, non-ground-truth (ngt) run list, then
    slides a fixed-duration window along each run's frames, computing the
    percent censored frames and mean FD/DV motion metrics within each
    window. Each output row corresponds to one window from one subID/run.
    Only subIDs that passed ground-truth requirements and have at least
    one qualifying window are included. Results are written to
    `window_list_1.txt` in `rep_sel_subID_window_dir`.

    Args:
        main_config_file (str): Path to the main pipeline configuration
            JSON file.
        pl_config_file (str): Path to the pipeline-specific configuration
            JSON file.
        preproc_type (str): Preprocessing type used to resolve each subID's
            preprocessed directory (e.g. "mc").
        run_type (str): Run type used to resolve each subID's preprocessed
            directory (e.g. "gt").
    """

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rep_sel_subID_ngt_dir = pl_config['path']['rep_sel_subID_ngt_dir']
    rep_sel_subID_window_dir = pl_config['path']['rep_sel_subID_window_dir']
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    FD_th = pl_config['preproc']['FD_th']
    DV_th = pl_config['preproc']['DV_th']
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    window_time_sec = pl_config['select_subID']['window']['window_time_sec']
    idx_frame_start_ini = pl_config['select_subID']['window'][
        'idx_frame_start_ini']
    col_ngt_nframes_name = pl_config['select_subID']['window'][
        'col_ngt_nframes_name']

    # Read input files
    ngt_run_list = os.path.join(rep_sel_subID_ngt_dir, 'ngt_list_1.txt')
    df = pd.read_csv(ngt_run_list, dtype={'bld': str})

    # Create output dirs
    if not os.path.exists(rep_sel_subID_window_dir):
        os.makedirs(rep_sel_subID_window_dir)

    # Initialize output df
    df_window = pd.DataFrame()
    for subID_tuple, df_subID in df.groupby(['dataset', 'subID']):
        # initialize subID-level variables
        windowID = 1
        idx_frame_start = idx_frame_start_ini
        for idx, row in df_subID.iterrows():
            # extract run-level data
            dataset = row['dataset']
            subID = row['subID']
            bld = row['bld']
            nframes_expected = row[col_ngt_nframes_name]
            TR = row['TR']

            # assign in dirs
            preproc_dir = util.get_preproc_dir_path(main_config_file,
                                                    pl_config_file,
                                                    preproc_type,
                                                    dataset,
                                                    subID,
                                                    run_type=run_type)
            assert os.path.isdir(
                preproc_dir), f"ERROR: {preproc_dir} does NOT exist!"

            # assign in files
            in_outliers_file = os.path.join(
                preproc_dir, 'qc',
                f"{subID}_bld{bld}{FD_stem}{FD_th}{DV_stem}{DV_th}_motion_outliers.txt"
            )
            assert os.path.isfile(
                in_outliers_file), f"ERROR: {in_outliers_file} does NOT exist!"
            df_outliers = pd.read_csv(in_outliers_file,
                                      header=None,
                                      names=['frames'],
                                      dtype={'bld': str})

            # Generate stem
            pre_mc_stem = util.generate_pre_mc_stem_from_preproc_config(
                pl_config_file, preproc_dir)

            # Read FD_metric outliers file
            in_outliers_FD_file = os.path.join(
                preproc_dir, 'bold', 'mc',
                f"{subID}_bld{bld}{pre_mc_stem}_motion_outliers{FD_stem}")
            assert os.path.isfile(
                in_outliers_FD_file
            ), f"ERROR: {in_outliers_FD_file} does NOT exist!"
            df_outliers_FD = pd.read_csv(in_outliers_FD_file,
                                         header=None,
                                         names=[f'{FD_metric}'],
                                         dtype={'bld': str})

            # Read DV_metric outliers file
            in_outliers_DV_file = os.path.join(
                preproc_dir, 'bold', 'mc',
                f"{subID}_bld{bld}{pre_mc_stem}_motion_outliers{DV_stem}")
            assert os.path.isfile(
                in_outliers_DV_file
            ), f"ERROR: {in_outliers_DV_file} does NOT exist!"
            df_outliers_DV = pd.read_csv(in_outliers_DV_file,
                                         header=None,
                                         names=[f'{DV_metric}'],
                                         dtype={'bld': str})

            # Retrieve number of frames for run
            nframes = len(df_outliers)
            assert nframes == nframes_expected, f"ERROR: (nframes={nframes}) != (nframes_expected={nframes_expected})"

            # Compute number of TRs per sliding window
            TRs_per_window = math.floor(window_time_sec / TR)

            # Loop across run-level; slide window within run
            for idx_window_start in range(idx_frame_start, nframes,
                                          TRs_per_window):
                idx_window_end = idx_window_start + TRs_per_window - 1

                # If window's last frame exceeds total number of frames, exit run-level loop
                if idx_window_end > nframes - 1:
                    break
                # print(subID_tuple, bld, nframes, idx_window_start, idx_window_end)

                # Compute start & end indices of window
                row_window = row
                if windowID < 10:
                    row_window['windowID'] = f'win-0{windowID}'
                else:
                    row_window['windowID'] = f'win-{windowID}'
                row_window['window_idx_start'] = idx_window_start
                row_window['window_idx_end'] = idx_window_end

                # Compute % censored frames of window
                df_window_outliers = df_outliers.iloc[
                    idx_window_start:idx_window_end + 1]
                assert len(df_window_outliers
                           ) == TRs_per_window, "ERROR: number of rows in \
                    df_window_outliers NOT equal to number of TRs per window"

                nframes_window = len(df_window_outliers)
                nframes_censored = (df_window_outliers['frames'] == 0).sum()
                nframes_noncensored = (df_window_outliers['frames'] == 1).sum()
                assert (
                    nframes_censored + nframes_noncensored
                ) == nframes_window, f"ERROR: censored + non-censored frames != total number of frames within window."
                percent_censored = (nframes_censored / nframes_window) * 100
                row_window['window_percent_censored'] = round(
                    np.float64(percent_censored), 4)

                # Compute mean FD_metric of window
                df_window_outliers_FD = df_outliers_FD.iloc[
                    idx_window_start:idx_window_end + 1]
                assert len(df_window_outliers_FD
                           ) == TRs_per_window, "ERROR: number of rows in \
                    df_window_outliers_FD NOT equal to number of TRs per window"

                window_mean_FD = df_window_outliers_FD.mean()
                row_window[f'window_mean_{FD_metric}'] = round(
                    np.float64(window_mean_FD.iloc[0]), 4)

                # Compute mean DV_metric of window
                df_window_outliers_DV = df_outliers_DV.iloc[
                    idx_window_start:idx_window_end + 1]
                assert len(df_window_outliers_DV
                           ) == TRs_per_window, "ERROR: number of rows in \
                    df_window_outliers_DV NOT equal to number of TRs per window"

                window_mean_DV = df_window_outliers_DV.mean()
                row_window[f'window_mean_{DV_metric}'] = round(
                    np.float64(window_mean_DV.iloc[0]), 4)

                # Append to output df_window
                df_row_window = pd.DataFrame(row_window).T
                df_window = pd.concat([df_window, df_row_window],
                                      ignore_index=True)

                # Update window ID
                windowID += 1

            # Compute subID-level % censored frames distribution
            if not df_window.empty:
                df_window_subID = df_window[(df_window['dataset'] == dataset)
                                            & (df_window['subID'] == subID)]

                # Check windowID uniqueness
                assert df_window_subID[
                    'windowID'].is_unique, f'ERROR: for {dataset}, {subID}, {bld}, windowIDs NOT unique!'

    # Save outputs
    df_window.to_csv(os.path.join(rep_sel_subID_window_dir,
                                  f'window_list_1.txt'),
                     index=False)


def plot_window_distribution(subID_list, windowID_list,
                             rep_sel_subID_window_dir):
    """Plot the distribution of window-level percent-censored-frames per subID.

    For every subID in `subID_list`, filters the window-level records in
    `windowID_list` down to that subID's windows, then plots and saves a
    histogram of `window_percent_censored` values across those windows.
    Each plot is written as a PNG to a `window_plot` subdirectory of
    `rep_sel_subID_window_dir`. If the `window_percent_censored` column
    is absent for a subID, a warning is printed and that subID is
    skipped.

    Args:
        subID_list (str): Path to a CSV file listing the subIDs to plot.
        windowID_list (str): Path to a CSV file (e.g. `window_list_1.txt`)
            containing per-window records, including
            `window_percent_censored`.
        rep_sel_subID_window_dir (str): Directory under which the
            `window_plot` subdirectory of output PNGs is created.
    """

    subID_df = pd.read_csv(subID_list, dtype={'bld': str})
    windowID_df = pd.read_csv(windowID_list, dtype={'bld': str})

    # Loop through each row in subID_df
    for tup in subID_df.iterrows():
        dataset = tup[1]['dataset']
        subID = tup[1]['subID']

        # Filter the windowID_df for matching dataset and subID
        windowID_subID_df = windowID_df[(windowID_df['dataset'] == dataset)
                                        & (windowID_df['subID'] == subID)]

        # Check if the 'window_percent_censored' column exists
        if 'window_percent_censored' in windowID_subID_df.columns:
            censored_values = windowID_subID_df['window_percent_censored']

            # Plot the distribution
            plt.figure(figsize=(10, 6))
            plt.hist(censored_values,
                     bins=range(0, 101, 10),
                     color='skyblue',
                     edgecolor='black',
                     alpha=0.7)
            plt.title(
                f"Distribution of {dataset} {subID}'s windows' mean % censored frames",
                fontsize=14)
            plt.xlabel("Mean % censored frames (within window)", fontsize=12)
            plt.ylabel("Frequency", fontsize=12)
            plt.xticks(range(0, 101, 10))

            # Set y-axis ticks to start from 0 and go in increments of 10
            plt.yticks(range(0,
                             int(max(plt.ylim())) + 10,
                             10))  # Calculate max y-limit and adjust y-ticks
            plt.grid(axis='y', alpha=0.75)

            # Define the output paths
            output_dir = os.path.join(rep_sel_subID_window_dir, "window_plot")
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, f"{dataset}_{subID}.png")

            # Save the plot
            plt.savefig(output_file)

            # Close the plot to avoid memory issues with many plots
            plt.close()

        else:
            print(
                f"Warning: 'window_percent_censored' column not found for {dataset} {subID}"
            )


def submit_select_subID_window_job(timestamp, jobprefix, main_config_file,
                                   pl_config_file, rep_sel_subID_ngt_dir,
                                   rep_sel_subID_window_dir, std_dir):
    """Submit a job to generate the window list and plot its distribution.

    Builds a Python one-liner that calls `generate_window_list` followed
    by `plot_window_distribution` with the appropriate preprocessing/run
    type ("mc"/"gt", mirroring `main`), then submits it to the job
    scheduler via `util.submit_job`.

    Args:
        timestamp (str): Timestamp string identifying this job submission
            (accepted for interface consistency with other job-submission
            functions; not otherwise used in this function's body).
        jobprefix (str): Prefix used to name the submitted job.
        main_config_file (str): Path to the main pipeline configuration
            JSON file.
        pl_config_file (str): Path to the pipeline-specific configuration
            JSON file.
        rep_sel_subID_ngt_dir (str): Directory containing the ngt subID
            list (`subID_list_1.txt`) used as input for plotting.
        rep_sel_subID_window_dir (str): Directory where the window list is
            written and the job's std/err logs are read from.
        std_dir (str): Directory to create for the job's stdout/stderr log
            files.
    """
    # Job directories and paths
    os.makedirs(std_dir, exist_ok=True)

    # Define types (mirrors `main`)
    preproc_type = "mc"
    run_type = "gt"

    # Define input/output files (mirrors `main`)
    subID_list = os.path.join(rep_sel_subID_ngt_dir, 'subID_list_1.txt')
    windowID_list = os.path.join(rep_sel_subID_window_dir, 'window_list_1.txt')

    # Generate job-specific variables
    jobname = f"{jobprefix}"
    stderr_file = os.path.join(std_dir, f"{jobname}.err")
    stdout_file = os.path.join(std_dir, f"{jobname}.out")

    # Construct command
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; sys.path.append(\\\"{TAN2026_TMSMOTION_DIR}/select_subID/window\\\"); "
        f"import CBIG_TMSmotion_select_subID_window as sel_subID_window; "
        f"sel_subID_window.generate_window_list(\\\"{main_config_file}\\\", \
            \\\"{pl_config_file}\\\", \\\"{preproc_type}\\\", \\\"{run_type}\\\"); "
        f"sel_subID_window.plot_window_distribution(\\\"{subID_list}\\\", \
            \\\"{windowID_list}\\\", \\\"{rep_sel_subID_window_dir}\\\")'; "
        f"conda deactivate")

    # Submit job to headnode
    util.submit_job(walltime="24:00:00",
                    memory="8G",
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def check_select_subID_window_job(timestamp, try_count,
                                  rep_sel_subID_window_dir):
    """Check whether the select_subID_window job completed successfully.

    Verifies that the expected output file (`window_list_1.txt`) exists
    in `rep_sel_subID_window_dir`. If it is missing, appends a failure
    message to a per-try log file (`fail_msg_try_<try_count>.txt`) in a
    timestamped subdirectory.

    Args:
        timestamp (str): Timestamp string used to name the output
            subdirectory where failure messages are logged.
        try_count (int): Current retry attempt count, used to name the
            failure message log file.
        rep_sel_subID_window_dir (str): Directory containing the job's
            expected output file(s).

    Returns:
        bool: True if all expected output files exist, False otherwise.
    """
    # Assume job successful
    all_job_success = True

    # Define output dir
    timestamp_dir = os.path.join(rep_sel_subID_window_dir, timestamp)
    os.makedirs(timestamp_dir, exist_ok=True)
    fail_msg = os.path.join(timestamp_dir, f"fail_msg_try_{try_count}.txt")

    # Check job output files for success
    filename2check = ["window_list_1.txt"]
    for filename in filename2check:
        exist = util.check_job_output(
            os.path.join(rep_sel_subID_window_dir, filename))
        if not exist:
            all_job_success = False
            # Handle failure message list
            if not os.path.exists(fail_msg):
                with open(fail_msg, 'w') as f:
                    f.write("error_message\n")
            with open(fail_msg, 'a') as f:
                f.write(f"{filename}\n")

    return all_job_success


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    """Generate the window list and plot its distribution for a pipeline run.

    Reads the main and pipeline-specific configuration files (pipeline
    config path given as the first command-line argument), then calls
    `generate_window_list` to build the sliding-window motion QC list
    for subIDs that passed ground-truth requirements, and
    `plot_window_distribution` to plot the resulting
    percent-censored-frames distributions.
    """
    # Read main config
    main_config_file = os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                                    'config', 'main_config.json')
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config_file = sys.argv[1]
    pl_config = util.read_json_file(pl_config_file)
    rep_sel_subID_ngt_dir = pl_config['path']['rep_sel_subID_ngt_dir']
    rep_sel_subID_window_dir = pl_config['path']['rep_sel_subID_window_dir']

    # Define types
    preproc_type = "mc"
    run_type = "gt"

    # Assign input files
    subID_list = os.path.join(rep_sel_subID_ngt_dir, 'subID_list_1.txt')
    windowID_list = os.path.join(rep_sel_subID_window_dir, 'window_list_1.txt')

    # Generate window list
    generate_window_list(main_config_file, pl_config_file, preproc_type,
                         run_type)

    # Plot distribution of subIDs who passed ground truth requirements
    plot_window_distribution(subID_list, windowID_list,
                             rep_sel_subID_window_dir)


if __name__ == "__main__":
    main()
