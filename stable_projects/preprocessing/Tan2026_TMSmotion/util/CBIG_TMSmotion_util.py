'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import json
import re
import subprocess
import time
import os
import shutil
import pandas as pd
import sys
import itertools

CBIG_CODE_DIR = os.environ.get('CBIG_CODE_DIR')


#############################################################################
# ===========================================================================
# GENERAL
# ===========================================================================
#############################################################################
def read_json_file(json_filename):
    """Read a JSON file and parse it into a Python dictionary.

    Args:
        json_filename (str): Path to the JSON file to read.

    Returns:
        dict: The parsed contents of the JSON file.
    """
    with open(json_filename, 'r') as json_file:
        data = json.load(json_file)
    return data


def drop_col(df, col2drop):
    """Drop the specified column(s) from a DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame.
        col2drop (list of str): Column names to drop. Example: ["nframes", "scan_time", "TR"].

    Returns:
        pd.DataFrame: The DataFrame with the specified columns removed (missing columns are ignored).
    """
    df = df.drop(columns=col2drop, errors='ignore')
    return df


def restart_dir(dir_path):
    """
    Deletes the directory if it exists and recreates it.

    Args:
        dir_path (str): Path to the directory to restart.
    """
    if os.path.islink(dir_path):
        os.unlink(dir_path)
    elif os.path.isdir(dir_path):
        shutil.rmtree(dir_path)

    os.makedirs(dir_path, exist_ok=True)


def write_config_dict2json(config, output_filename):
    """
    Function to write the configuration dictionary to a JSON file.

    Args:
        config (dict): The configuration dictionary to be saved.
        output_filename (str): The full path for output JSON.
    """
    # Create parent directories if they don't exist
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    # Write the dictionary to a JSON file
    with open(output_filename, 'w') as json_file:
        json.dump(config, json_file, indent=4)

    # Print success message
    print(
        f"JSON configuration file '{output_filename}' generated successfully.")


def save_df_to_csv(df, output_dir, filename):
    """
    Save a DataFrame to a CSV file in the specified directory.

    Args:
        df (pd.DataFrame): The DataFrame to save.
        output_dir (str): The directory where the file will be saved.
        filename (str): The name of the CSV file.

    Returns:
        str: The full path of the saved file.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, filename)
    df.to_csv(output_path, index=False)
    return output_path


def filter_dataset_subID_from_all_subID_df(df, dataset, subID):
    """
    Filters the input dataframe based on the provided dataset and subID.

    Args:
        df (pd.DataFrame): The dataframe to filter.
        dataset (str): The dataset value to filter by.
        subID (str): The subID value to filter by.

    Returns:
        pd.DataFrame: A filtered dataframe containing only rows with the specified dataset and subID.
    """
    filtered_df = df[(df['dataset'] == dataset) & (df['subID'] == subID)]

    return filtered_df


def read_filter_save_selected_subID_list(base_dir, filename, dataset, subID):
    """
    Reads a CSV file, filters the rows based on the dataset and subID,
    and saves the filtered DataFrame in a subdirectory based on dataset and subID.

    Args:
        base_dir (str): The base directory where the CSV file is located and
                        where the filtered output will be saved.
        filename (str): The name of the file (e.g., "gt_list_1.txt") to read
                        and save. This file is expected to be in CSV format.
        dataset (str): The dataset value to filter by (e.g., 'dataset1').
        subID (str): The subject ID value to filter by (e.g., 'sub01').

    Returns:
        str: The file path of the saved filtered CSV file, which includes the
             dataset and subID-specific subdirectory.
    """

    # Read the DataFrame
    input_file_path = os.path.join(base_dir, filename)
    df = pd.read_csv(input_file_path,
                     dtype={
                         'rm_run_th': str,
                         'bld': str,
                         'sim_bld': str,
                         'highest_sim_bld': str
                     })

    # Filter based on dataset and subID
    filtered_df = filter_dataset_subID_from_all_subID_df(df, dataset, subID)

    # Save filtered df
    base_filename, _ = os.path.splitext(filename)
    output_dir = os.path.join(base_dir, base_filename, dataset, subID)
    save_df_to_csv(filtered_df, output_dir, filename)

    # Construct output file path
    output_file_path = os.path.join(output_dir, filename)

    return output_file_path


def reorder_df_columns(df, desired_columns):
    """
    Reorders the columns of a DataFrame based on a desired list of column names.

    Args:
        df (pd.DataFrame): The DataFrame whose columns need to be reordered.
        desired_columns (list of str): A list of column names specifying the desired order.
            Columns not in this list will be appended at the end.

    Returns:
        pd.DataFrame: A new DataFrame with reordered columns.
    """
    # Separate desired columns and remaining columns
    remaining_columns = [
        col for col in df.columns if col not in desired_columns
    ]
    new_column_order = desired_columns + remaining_columns

    # Reorder the DataFrame
    return df[new_column_order]


def delete_file(file2del):
    """
    Deletes the specified file if it exists.

    Args:
        file2del (str): Path to the file to be deleted.
    """
    if os.path.exists(file2del):
        print(f"Deleting existing file: {file2del}")
        os.remove(file2del)


def filter_sim_runID(subID_df, num_sim_run_per_sim_ses):
    """
    Filters a DataFrame to retain rows where the sim_runID is within the allowed range.

    Args:
        subID_df (pd.DataFrame): The input DataFrame containing a 'sim_runID' column.
        num_sim_run_per_sim_ses (int): The maximum number of simulation runs per session to retain.

    Returns:
        pd.DataFrame: A filtered DataFrame with rows containing sim_runIDs within the allowed range.
    """
    # Extract the numeric part of 'sim_runID' and filter based on num_sim_run_per_sim_ses
    subID_df['sim_run_number'] = subID_df['sim_runID'].str.extract(
        r'run-(\d+)').astype(int)
    filtered_df = subID_df[subID_df['sim_run_number'] <=
                           num_sim_run_per_sim_ses].drop(
                               columns='sim_run_number')
    return filtered_df


#############################################################################
# ===========================================================================
# ENVIRONMENT VARIABLE GENERATION
# ===========================================================================
#############################################################################


def get_cbig_code_dir():
    """
    Retrieve the CBIG_CODE_DIR environment variable and ensure it is set.

    Returns:
        str: The value of the CBIG_CODE_DIR environment variable.

    Raises:
        EnvironmentError: If CBIG_CODE_DIR is not set.
    """
    cbig_code_dir = os.getenv('CBIG_CODE_DIR')
    if not cbig_code_dir:
        raise EnvironmentError(
            "CBIG_CODE_DIR environment variable is not set.")
    return cbig_code_dir


#############################################################################
# ===========================================================================
# STEM GENERATION
# ===========================================================================
#############################################################################
def generate_stem_from_config(pl_config_file, config_file):
    """Generate the BOLD filename stem based on a CBIG preproc config file.

    Parses the CBIG preproc config file line by line, appending a suffix to
    the stem for each recognized preprocessing step encountered (skip, slice
    timing correction, motion correction with FD/DV thresholds, nuisance
    regression, and bandpass filtering), building up the stem incrementally.

    Args:
        pl_config_file (str): Path to the pipeline configuration JSON file
            (used to look up FD_stem, DV_stem, and bandpass frequencies).
        config_file (str): Path to the CBIG preproc config file to parse.

    Returns:
        tuple: (final_stem, stem_list)
            - final_stem (str): The fully built stem after all recognized
              steps have been applied.
            - stem_list (list of str): The stem after each recognized step,
              in the order the steps were encountered.
    """
    # Read pipeline config
    pl_config = read_json_file(pl_config_file)
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    bandpass_low_f = pl_config['preproc']['bandpass_low_f']
    bandpass_high_f = pl_config['preproc']['bandpass_high_f']

    # Init final_stem
    final_stem = '_rest'
    stem_list = ['_rest']

    # Validate the config file exists
    if not os.path.exists(config_file):
        raise FileNotFoundError(f'ERROR: Config file not found: {config_file}')

    # Parse the config file
    try:
        with open(config_file, 'r') as f:
            for line in f:
                # Check for CBIG_preproc_skip
                if "CBIG_preproc_skip" in line:
                    skip_value = line.split("-skip")[-1].strip().split()[0]
                    final_stem += f"_skip{skip_value}"
                    stem_list.append(final_stem)
                # Check for CBIG_preproc_fslslicetimer
                if "CBIG_preproc_fslslicetimer" in line:
                    final_stem += "_stc"
                    stem_list.append(final_stem)
                # Check for CBIG_preproc_fslmcflirt_outliers
                if "CBIG_preproc_fslmcflirt_outliers" in line:
                    final_stem += "_mc"
                    stem_list.append(final_stem)

                    # Extract FD_th and DV_th from config
                    FD_th_match = re.search(r'-FD_th\s+([\d\.]+)', line)
                    DV_th_match = re.search(r'-DV_th\s+([\d\.]+)', line)
                    if FD_th_match:
                        FD_th = FD_th_match.group(1)
                    if DV_th_match:
                        DV_th = DV_th_match.group(1)
                # Check for CBIG_preproc_regress
                if "CBIG_preproc_regress" in line:
                    final_stem += f"_residc_interp{FD_stem}{FD_th}{DV_stem}{DV_th}"
                    stem_list.append(final_stem)
                # Check for CBIG_preproc_fslmcflirt_outliers
                if "CBIG_preproc_bandpass" in line:
                    final_stem += f"_bp_{bandpass_low_f}_{bandpass_high_f}"
                    stem_list.append(final_stem)
    except IOError as e:
        raise IOError(
            f"ERROR: Failed to read config file: {config_file}") from e

    return final_stem, stem_list


def generate_pre_mc_stem_from_precursor_stem(stem):
    """
    Generate the pre-mc stem by removing everything after '_mc' in the given stem string.

    Args:
        stem (str): The stem string to be processed.

    Returns:
        str: The modified stem string without '_mc' and subsequent parts.

    Raises:
        ValueError: If '_mc' is not found in the stem.
    """
    # Find the position of '_mc' in the stem
    if "_mc" in stem:
        # Cut off everything after '_mc' (including '_mc')
        pre_mc_stem = stem.split("_mc")[0]
        return pre_mc_stem
    else:
        raise ValueError(f"'_mc' could not be found in {stem} stem.")


def generate_pre_mc_stem_from_preproc_config(pl_config_file, preproc_dir):
    """
    Generate the preprocessing motion correction stem based on a configuration file.

    Args:
        pl_config_file (str): Path to the pipeline configuration JSON file.
        preproc_dir (str): Path to the preprocessing directory containing the configuration file.

    Returns:
        str: The motion correction stem extracted from the configuration file.
    """
    # Read pipeline config
    read_json_file(pl_config_file)

    # Define the path to the configuration file
    config_file = os.path.join(preproc_dir, 'logs',
                               'CBIG_preproc_fMRI_preprocess.config')

    # Generate the resume stem from the configuration file
    final_stem, _ = generate_stem_from_config(pl_config_file, config_file)

    # Generate the motion correction stem from the resume stem
    pre_mc_stem = generate_pre_mc_stem_from_precursor_stem(final_stem)

    return pre_mc_stem


def append_suffix_to_stem(stem, suffix):
    """
    Appends a suffix to a given string stem.

    Args:
        stem (str): The original string to which the suffix will be appended.
        suffix (str): The suffix to append to the stem.

    Returns:
        str: The updated string with the suffix appended.

    Raises:
        ValueError: If stem is not a string.
    """
    if not isinstance(stem, str):
        raise ValueError("The 'stem' argument must be a string.")
    return stem + suffix


def generate_pre_FD_and_post_DV_stem(BOLD_stem, FD_stem, DV_stem):
    """Generate pre_FD_stem and post_DV_stem from a BOLD stem string.

    Args:
        BOLD_stem (str): The full BOLD stem string to split.
        FD_stem (str): The FD (framewise displacement) marker substring
            to locate within BOLD_stem.
        DV_stem (str): The DV (DVARS) marker substring to locate within
            BOLD_stem.

    Returns:
        tuple: (pre_FD_stem, post_DV_stem)
            - pre_FD_stem (str): The substring before FD_stem.
              Example: "_rest_skip10_stc_mc_residc_interp"
            - post_DV_stem (str): The substring after the DV metric value and includes additional parts.
              Example: "_bp_0.009_0.08"

    Raises:
        ValueError: If BOLD_stem does not contain both FD_stem and DV_stem.
    """
    # Find the indices of FD_stem and DV_stem
    fd_index = BOLD_stem.find(FD_stem)
    dv_index = BOLD_stem.find(DV_stem)

    if fd_index == -1 or dv_index == -1:
        raise ValueError("Input string must contain both FD_stem and DV_stem.")

    # Extract pre_FD_stem
    pre_FD_stem = BOLD_stem[:fd_index]

    # Extract post_DV_stem
    post_DV_start = BOLD_stem.find(
        '_', dv_index + len(DV_stem))  # Find the underscore after DV metric
    post_DV_stem = BOLD_stem[post_DV_start:]

    return pre_FD_stem, post_DV_stem


#############################################################################
# ===========================================================================
# PATH GENERATION
# ===========================================================================
#############################################################################
def get_reconall_base_dir_path(reconall_dataset, reconall_subID,
                               reconall_sesID, data_dir):
    """
    Generate the full path for the reconall base directory.

    Args:
        reconall_dataset (str): The dataset name for reconall.
        reconall_subID (str): The subject ID for reconall.
        reconall_sesID (str): The session ID for reconall.
        data_dir (str): The base directory where the data is stored.

    Returns:
        str: The full path to the reconall base directory.
    """
    # Generate the full path using os.path.join
    base_dir_path = os.path.join(data_dir, reconall_dataset, 'preprocessed',
                                 'reconall', reconall_subID, reconall_sesID)

    return base_dir_path


def get_betsurf_dir_path(dataset, subID, sesID, runID, data_dir):
    """
    Generate the full path for the betsurf base directory.

    Args:
        dataset (str): The dataset name for betsurf.
        subID (str): The subject ID for betsurf.
        sesID (str): The session ID for betsurf.
        runID (str): The run ID for betsurf.
        data_dir (str): The base directory where the data is stored.

    Returns:
        str: The full path to the betsurf base directory.
    """
    # Generate the full path using os.path.join
    dir_path = os.path.join(data_dir, dataset, 'preprocessed', 'betsurf',
                            subID, sesID, runID)

    return dir_path


def get_preproc_dir_path(main_config_file,
                         pl_config_file,
                         preproc_type,
                         dataset,
                         subID,
                         sesID=None,
                         runID=None,
                         run_type=None):
    """Generate the full path to the preprocessing directory.

    The resulting path layout depends on preproc_type: 'before_mc' paths
    are keyed by sesID/runID, while 'mc' and 'after_mc' paths are keyed by
    run_type instead.

    Args:
        main_config_file (str): Path to the main configuration JSON file.
        pl_config_file (str): Path to the pipeline configuration JSON file.
        preproc_type (str): One of 'before_mc', 'mc', or 'after_mc'.
        dataset (str): The dataset identifier.
        subID (str): The subject ID.
        sesID (str, optional): The session ID. Required when preproc_type
            is 'before_mc'. Defaults to None.
        runID (str, optional): The run ID. Required when preproc_type is
            'before_mc'. Defaults to None.
        run_type (str, optional): The run type (e.g., 'ngt' or 'gt').
            Required when preproc_type is 'mc' or 'after_mc'. Defaults to None.

    Returns:
        str: The full path to the preprocessing directory.

    Raises:
        ValueError: If preproc_type is not one of 'before_mc', 'mc', or 'after_mc'.
    """
    # Read main config
    main_config = read_json_file(main_config_file)
    data_dir = main_config['path']['data_dir']

    # Read pipeline config
    pl_config = read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])

    # Define preproc dir path
    if preproc_type == "before_mc":
        preproc_dir = os.path.join(data_dir, dataset, 'preprocessed', 'rest',
                                   pl_name, preproc_type, subID, sesID, runID)
    elif preproc_type in {"mc", "after_mc"}:
        preproc_dir = os.path.join(data_dir, dataset, 'preprocessed', 'rest',
                                   pl_name, f"{preproc_type}_{run_type}",
                                   subID)
    else:
        raise ValueError(
            f"Invalid preproc_type '{preproc_type}'. Accepted values are 'before_mc', "
            "'mc', or 'after_mc'.")

    return preproc_dir


def get_sim_ses_preproc_dir_real_path(main_config_file, pl_config_file,
                                      dataset, subID, sim_run_bin, sim_sesID):
    """
    Generate full path to base directory storing actual simulation run
    instances NIFTIs.

    Args:
        main_config_file (str): Path to the main configuration JSON file.
        pl_config_file (str): Path to the pipeline configuration JSON file.
        dataset (str): The dataset identifier.
        subID (str): The subject ID.
        sim_run_bin (str): The simulation run bin identifier.
        sim_sesID (str): The simulation session ID.

    Returns:
        str: The full path to the preprocessed base directory.
    """
    # Read main config
    main_config = read_json_file(main_config_file)
    rep_dir = main_config['path']['rep_dir']

    # Read pipeline config
    pl_config = read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])

    # Generate the full path using os.path.join
    sim_ses_preproc_real_dir = os.path.join(rep_dir, pl_name, 'output',
                                            'real_sim', dataset, subID,
                                            sim_run_bin, sim_sesID)

    return sim_ses_preproc_real_dir


def get_sim_ses_preproc_dir_symlink_path(main_config_file,
                                         pl_config_file,
                                         dataset,
                                         subID,
                                         rm_run_th=None,
                                         combo=None,
                                         sim_sesID=None):
    """
    Generate full path to base directory storing symlink simulation run
    instances NIFTIs.

    If rm_run_th, combo, and sim_sesID are all provided, the path is
    narrowed down to that specific combination; otherwise the base
    per-subject symlink directory is returned.

    Args:
        main_config_file (str): Path to the main configuration JSON file.
        pl_config_file (str): Path to the pipeline configuration JSON file.
        dataset (str): The dataset identifier.
        subID (str): The subject ID.
        rm_run_th (str, optional): Removed-run threshold identifier. Defaults to None.
        combo (str, optional): The FD/DV combination identifier. Defaults to None.
        sim_sesID (str, optional): The simulation session ID. Defaults to None.

    Returns:
        str: The full path to the preprocessed base directory.
    """
    # Read main config
    main_config = read_json_file(main_config_file)
    rep_dir = main_config['path']['rep_dir']

    # Read pipeline config
    pl_config = read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])

    # Generate the full path using os.path.join
    if rm_run_th and combo and sim_sesID:
        sim_ses_preproc_symlink_dir = os.path.join(rep_dir, pl_name, 'output',
                                                   'symlink_sim', dataset,
                                                   subID, rm_run_th, combo,
                                                   sim_sesID)
    else:
        sim_ses_preproc_symlink_dir = os.path.join(rep_dir, pl_name, 'output',
                                                   'symlink_sim', dataset,
                                                   subID)

    return sim_ses_preproc_symlink_dir


def get_split_bld_preproc_dir_path(main_config_file, pl_config_file, dataset,
                                   subID):
    """
    Generate the full path for the ground truth
    splitted preprocessed directory.

    Args:
        main_config_file (str): Path to the main configuration JSON file.
        pl_config_file (str): Path to the pipeline configuration JSON file.
        dataset (str): The dataset identifier.
        subID (str): The subject ID.

    Returns:
        str: The full path to the preprocessed base directory.
    """

    # Read main config
    main_config = read_json_file(main_config_file)
    rep_dir = main_config['path']['rep_dir']

    # Read pipeline config
    pl_config = read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])

    # Generate the full path using os.path.join
    split_bld_preproc_dir = os.path.join(rep_dir, pl_name, 'output',
                                         'gt_split', dataset, subID)

    return split_bld_preproc_dir


def get_config_file_dir_path(main_config_file,
                             pl_config_file,
                             timestamp,
                             preproc_type,
                             dataset,
                             subID,
                             sesID=None,
                             runID=None,
                             run_type=None):
    """Generate the full path to the preproc config.txt file.

    The resulting path layout depends on preproc_type: 'before_mc' paths
    are keyed by sesID/runID, while 'mc' and 'after_mc' paths are keyed by
    run_type instead.

    Args:
        main_config_file (str): Path to the main configuration JSON file.
        pl_config_file (str): Path to the pipeline configuration JSON file
            (used to look up config_dir).
        timestamp (str): Timestamp identifier accepted by this function but
            not currently used in its body.
        preproc_type (str): One of 'before_mc', 'mc', or 'after_mc'.
        dataset (str): The dataset identifier.
        subID (str): The subject ID.
        sesID (str, optional): The session ID. Required when preproc_type
            is 'before_mc'. Defaults to None.
        runID (str, optional): The run ID. Required when preproc_type is
            'before_mc'. Defaults to None.
        run_type (str, optional): The run type (e.g., 'ngt' or 'gt').
            Required when preproc_type is 'mc' or 'after_mc'. Defaults to None.

    Returns:
        str: The full path to the config.txt file.

    Raises:
        ValueError: If preproc_type is not one of 'before_mc', 'mc', or 'after_mc'.
    """
    # Read main config
    read_json_file(main_config_file)

    # Read pipeline config
    pl_config = read_json_file(pl_config_file)
    config_dir = pl_config['path']['config_dir']

    if preproc_type == "before_mc":
        assert sesID is not None and runID is not None, (
            f"ERROR: for {preproc_type}, both 'sesID' and/or 'runID' cannot be None."
        )
        config_file = os.path.join(config_dir, preproc_type, dataset, subID,
                                   sesID, runID, "config.txt")
    elif preproc_type in {"mc", "after_mc"}:
        assert run_type is not None, f"ERROR: for {preproc_type} 'run_type' cannot be None."
        config_file = os.path.join(config_dir, f"{preproc_type}_{run_type}",
                                   dataset, subID, "config.txt")
    else:
        raise ValueError(
            f"Invalid preproc_type '{preproc_type}'. Accepted values are 'before_mc', "
            "'mc', or 'after_mc'.")

    return config_file


def get_TMS_dir_path(run_type,
                     main_config_file,
                     dataset,
                     pl_name,
                     subID,
                     rm_run_th=None,
                     combo=None,
                     sim_sesID=None,
                     FD_metric=None,
                     vary_FD_th=None,
                     DV_metric=None,
                     vary_DV_th=None,
                     TMS_type='Kong2025_dep'):
    """Generate the TMS_Kong2025 directory path.

    For run_type 'ngt' (non-ground-truth), the path is extended with the
    removed-run threshold, combo, simulation session, and FD/DV threshold
    subdirectories. For run_type 'gt' (ground-truth), a simpler 'gt' path
    is returned.

    Args:
        run_type (str): Either 'ngt' or 'gt'.
        main_config_file (str): Path to the main configuration JSON file.
        dataset (str): The dataset identifier.
        pl_name (str): The pipeline name.
        subID (str): The subject ID.
        rm_run_th (str, optional): Removed-run threshold identifier.
            Required when run_type is 'ngt'. Defaults to None.
        combo (str, optional): The FD/DV combination identifier. Required
            when run_type is 'ngt'. Defaults to None.
        sim_sesID (str, optional): The simulation session ID. Required
            when run_type is 'ngt'. Defaults to None.
        FD_metric (str, optional): The FD metric name. Required when
            run_type is 'ngt'. Defaults to None.
        vary_FD_th (str, optional): The FD threshold value being varied.
            Required when run_type is 'ngt'. Defaults to None.
        DV_metric (str, optional): The DV metric name. Required when
            run_type is 'ngt'. Defaults to None.
        vary_DV_th (str, optional): The DV threshold value being varied.
            Required when run_type is 'ngt'. Defaults to None.
        TMS_type (str, optional): The TMS output subdirectory name.
            Defaults to 'Kong2025_dep'.

    Returns:
        str: The full path to the TMS_Kong2025 directory.

    Raises:
        ValueError: If run_type is not 'ngt' or 'gt'.
    """
    # Read main config
    main_config = read_json_file(main_config_file)
    rep_dir = main_config['path']['rep_dir']

    # Core path for non-ground truth & ground truth cases
    base_dir = os.path.join(rep_dir, pl_name, 'output', f"TMS_{TMS_type}",
                            dataset, subID)

    # Append depending on non-ground truth or ground truth
    if run_type == "ngt":
        return os.path.join(
            base_dir, "ngt", rm_run_th, combo, sim_sesID,
            f"{FD_metric}{vary_FD_th}_{DV_metric}{vary_DV_th}")
    elif run_type == "gt":
        return os.path.join(base_dir, "gt")
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted values are 'ngt' or 'gt'."
        )


def get_compare_Kong2025_output_dir_path(main_config_file,
                                         pl_config_file,
                                         output_type,
                                         dataset,
                                         subID,
                                         sub_output_type=None,
                                         rm_run_th=None,
                                         combo=None,
                                         sim_sesID=None,
                                         FD_metric=None,
                                         vary_FD_th=None,
                                         DV_metric=None,
                                         vary_DV_th=None,
                                         TMS_type='Kong2025_dep'):
    """Construct the compare-Kong2025 output directory path.

    If sub_output_type is given, it is included in the path. If
    rm_run_th, combo, sim_sesID, FD_metric, vary_FD_th, DV_metric, and
    vary_DV_th are all provided, the path is extended with the
    non-ground-truth ('ngt') subdirectories; otherwise a 'gt'
    (ground-truth) subdirectory is used.

    Args:
        main_config_file (str): Path to the main configuration JSON file
            (currently read but unused in this function's body).
        pl_config_file (str): Path to the pipeline configuration JSON file
            (used to look up rep_compare_Kong2025_dep_dir).
        output_type (str): The output type subdirectory name.
        dataset (str): The dataset identifier.
        subID (str): The subject ID.
        sub_output_type (str, optional): A sub-output-type subdirectory
            name. Defaults to None.
        rm_run_th (str, optional): Removed-run threshold identifier. Defaults to None.
        combo (str, optional): The FD/DV combination identifier. Defaults to None.
        sim_sesID (str, optional): The simulation session ID. Defaults to None.
        FD_metric (str, optional): The FD metric name. Defaults to None.
        vary_FD_th (str, optional): The FD threshold value being varied. Defaults to None.
        DV_metric (str, optional): The DV metric name. Defaults to None.
        vary_DV_th (str, optional): The DV threshold value being varied. Defaults to None.
        TMS_type (str, optional): The TMS type, used to select the
            config-file path key. Defaults to 'Kong2025_dep'.

    Returns:
        str: The constructed output directory path.

    Raises:
        ValueError: If TMS_type is not 'Kong2025_dep'.
    """
    # Read main config
    read_json_file(main_config_file)

    # Read pipeline config
    pl_config = read_json_file(pl_config_file)
    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )

    # If sub_output_type is given, include it in the path
    if sub_output_type is not None:
        output_dir = os.path.join(rep_compare_Kong2025_dir, output_type,
                                  sub_output_type, dataset, subID)
    else:
        output_dir = os.path.join(rep_compare_Kong2025_dir, output_type,
                                  dataset, subID)

    # If rm_run_th, combo, sim_sesID, FD_metric, vary_FD_th, DV_metric, vary_DV_th are all provided, add to path
    if all(x is not None for x in [
            rm_run_th, combo, sim_sesID, FD_metric, vary_FD_th, DV_metric,
            vary_DV_th
    ]):
        output_dir = os.path.join(
            output_dir, "ngt", rm_run_th, combo, sim_sesID,
            f"{FD_metric}{vary_FD_th}_{DV_metric}{vary_DV_th}")
    else:
        output_dir = os.path.join(output_dir, "gt")

    return output_dir


#############################################################################
# ===========================================================================
# RECONALL IDENTITY EXTRACTION
# ===========================================================================
#############################################################################
def extract_reconall_map(reconall_map, dataset, subID):
    """
    Extracts the reconall mapping for the given dataset and subID.
    Rationale for function is because certain dataset's subID
    missing T1 images (i.e., Newbold, sub-cast1). However, because
    participant identity is the same as another dataset's subID
    (i.e. MSC, sub-MSC02), so can still use other dataset's T1 for
    same participant.

    Args:
        reconall_map (list): List of dictionaries containing the reconall map.
        dataset (str): The dataset to search for.
        subID (str): The subID to search for.

    Returns:
        tuple: A tuple (reconall_dataset, reconall_subID) if found, else None.
    """
    # Iterate through the reconall_map list
    for entry in reconall_map:
        # Check if both the dataset and subID match
        if entry['dataset'] == dataset and entry['subID'] == subID:
            # Return the mapped dataset and subID from the "map2" field
            reconall_dataset = entry['map2']['dataset']
            reconall_subID = entry['map2']['subID']
            return reconall_dataset, reconall_subID

    # If no match found, return None
    return dataset, subID


def get_first_reconall_sesID_runID(reconall_dataset, reconall_subID,
                                   first_N_T1_list):
    """
    Returns the sesID and runID of the first row within a group based on reconall_dataset and reconall_subID.
    If the group is not found or if sesID and runID are missing, an error will be raised.

    Args:
        reconall_dataset (str): The reconall dataset name.
        reconall_subID (str): The reconall subject ID.
        first_N_T1_list (pd.DataFrame): DataFrame containing the first_N_T1_list.

    Returns:
        tuple: (sesID, runID) of the first row for the specific reconall_dataset and reconall_subID group.

    Raises:
        ValueError: If the reconall_dataset and reconall_subID group is not found or missing sesID and runID.
    """
    # Group by reconall_dataset, reconall_subID
    df = pd.read_csv(first_N_T1_list)
    grouped = df.groupby(['dataset', 'subID'])

    # Check if the group exists
    if (reconall_dataset, reconall_subID) in grouped.groups:
        # Get the first row in that group
        group = grouped.get_group((reconall_dataset, reconall_subID))

        # Ensure sesID and runID are present
        if 'sesID' in group.columns and 'runID' in group.columns:
            # Return the sesID and runID of the first row
            sesID = group.iloc[0]['sesID']
            runID = group.iloc[0]['runID']
            return sesID, runID
        else:
            raise ValueError(
                f"Missing 'sesID' or 'runID' for {reconall_dataset} & {reconall_subID} in `group`."
            )

    # Raise an error if the group does not exist
    raise ValueError(
        f"Missing {reconall_dataset} & {reconall_subID} in first_N_T1_list.")


#############################################################################
# ===========================================================================
# RECONALL IDENTITY EXTRACTION
# ===========================================================================
#############################################################################
def generate_FD_DV_combination(FD_metric_value, DV_metric_value, mode):
    """
    Generate combinations of FD and DV values based on the specified mode.

    Args:
        FD_metric_value (list): A list of FD values to vary.
        DV_metric_value (list): A list of DV values to vary.
        mode (str): Mode of combination.
            "grid" for full Cartesian product,
            "slide" for slide combinations (1:1 pairing).

    Returns:
        list of tuples: Each tuple contains an (FD, DV) pair.

    Raises:
        ValueError: If mode is not "grid" or "slide".
    """
    if mode == "grid":
        # Full Cartesian product
        FD_DV_combination_list = list(
            itertools.product(FD_metric_value, DV_metric_value))
    elif mode == "slide":
        # Sliding combinations
        FD_DV_combination_list = [
            (FD_metric_value[i], DV_metric_value[i])
            for i in range(min(len(FD_metric_value), len(DV_metric_value)))
        ]
    else:
        raise ValueError("Invalid mode. Choose 'grid' or 'slide'.")

    return FD_DV_combination_list


#############################################################################
# ===========================================================================
# JOB SUBMISSION
# ===========================================================================
#############################################################################
def check_job_output(file):
    """
    Check if a file exists.

    Args:
        file (str): Path to the file to be checked.

    Returns:
        bool: True if the file exists, False otherwise.
    """
    if os.path.exists(file):
        return True
    else:
        return False


def check_job_text(file, text):
    """
    Check if the file exists and contains the given text.

    Args:
        file (str): Path to the file to be checked.
        text (str): The text to search for in the file.

    Returns:
        tuple: (bool, str)
               - True if the text is found in the file, False otherwise.
               - Message explaining the result.
    """
    if os.path.exists(file):
        # Check if the text is present in the file
        with open(file, 'r') as f:
            content = f.read()
            if text in content:
                return True, f"Text '{text}' found in {file}."
            else:
                return False, f"{file} exists, but '{text}' not found in file."
    else:
        return False, f"File '{file}' does not exist."


def check_job_status(jobprefix, username):
    """Check the status of jobs in the PBS queue with a specific jobprefix.

    Queries `qstat` for the given user and inspects the job name column for
    jobprefix, treating jobs marked ' R ' (running) or ' Q ' (queued) as
    still in progress.

    Args:
        jobprefix (str): Substring to match against job names.
        username (str): Username to query in qstat.

    Returns:
        bool: True if all jobs with the given prefix are finished (neither
            running nor queued), or if the qstat query fails. False if any
            matching job is still running or queued.
    """
    # Get the list of jobs from qstat
    qstat_command = f"qstat -u {username}"
    try:
        # Run the qstat command and capture the output
        result = subprocess.run(qstat_command,
                                shell=True,
                                capture_output=True,
                                text=True,
                                check=True)
        jobs_info = result.stdout

        # Check if there are any jobs with the jobprefix in running or queued state
        jobs_running_or_queued = [
            line for line in jobs_info.splitlines()
            if jobprefix in line and (' R ' in line or ' Q ' in line)
        ]
        print(f'{len(jobs_running_or_queued)} job(s) left...')
        # If no running or queued jobs are found, then all jobs have finished
        if len(jobs_running_or_queued) == 0:
            return True
        else:
            return False

    except subprocess.CalledProcessError as e:
        print(f"Error checking jobs status: {e}")
        return False


def wait_for_job(jobprefix, username):
    """Block until all jobs with the given jobprefix have finished.

    Polls check_job_status every 60 seconds until it reports all jobs
    with jobprefix are no longer running or queued.

    Args:
        jobprefix (str): Substring to match against job names.
        username (str): Username to query in qstat.
    """
    print(f"Waiting for {jobprefix} jobs to finish...")

    while not check_job_status(jobprefix, username):
        print(f"{jobprefix} jobs still running...")
        time.sleep(60)

    print(f"{jobprefix} jobs finished.")


def submit_job(walltime, memory, ncpus, jobname, stderr_file, stdout_file,
               cmd_script):
    """
    Submit a job using modified pbsubmit script.
    Modification is the removal of `ssh headnode` commands.
    Reason is that each ssh into headnode will generate a log file.
    If generate too many will fill up the log folder too fast.

    Args:
        walltime (str): Walltime for the job (e.g., "168:00:00").
        memory (str): Memory allocation (e.g., "20G").
        ncpus (int): Number of CPUs.
        jobname (str): Job name.
        stderr_file (str): Path to the stderr log file.
        stdout_file (str): Path to the stdout log file.
        cmd_script (str): Command to execute.
    """
    # Create parent dir(s) for std files
    os.makedirs(os.path.dirname(stderr_file), exist_ok=True)
    os.makedirs(os.path.dirname(stdout_file), exist_ok=True)
    # Prepare submit command
    CBIG_pbsubmit = os.path.join(CBIG_CODE_DIR, "setup", "CBIG_pbsubmit")
    qsub_command = (
        f"{CBIG_pbsubmit} -cmd \"{cmd_script}\" -walltime {walltime} -mem {memory} -ncpus {ncpus} "
        f"-name {jobname} -joberr {stderr_file} -jobout {stdout_file}")
    # Submit job
    subprocess.run(qsub_command, shell=True, check=True)
    # print(f"Submitted job {jobname}.")
    # print(qsub_command)
    # time.sleep(3)


def submit_TMS_Kong2025_job(walltime, memory, ncpus, jobname, stderr_file,
                            stdout_file, cmd_script, max_jobs, username,
                            jobprefix):
    """
    Submit a TMS_Kong2025 job using pbsubmit script.
    Rationale is because each job generates temporary ~20 GB nii file.
    Due to storage resource constraints, need to curb number of
    parallel running TMS_Kong2025 jobs.

    Args:
        walltime (str): Walltime for the job (e.g., "168:00:00").
        memory (str): Memory allocation (e.g., "20G").
        ncpus (int): Number of CPUs.
        jobname (str): Job name.
        stderr_file (str): Path to the stderr log file.
        stdout_file (str): Path to the stdout log file.
        cmd_script (str): Command to execute.
        max_jobs (int): Maximum number of concurrent jobs with jobprefix
            allowed for username before waiting to submit.
        username (str): Username to query in qstat.
        jobprefix (str): Substring to match against job names when
            counting currently running/queued jobs.
    """
    # Create parent dir(s) for std files
    os.makedirs(os.path.dirname(stderr_file), exist_ok=True)
    os.makedirs(os.path.dirname(stdout_file), exist_ok=True)

    # Wait for number of jobs with jobprefix & username to decline
    # below threshold before submitting new jobs
    while get_number_of_job(username, jobprefix) >= max_jobs:
        print(
            f"Waiting: Current {jobprefix} jobs for {username} >= {max_jobs}. Retrying in 10 seconds..."
        )
        time.sleep(10)

    # Prepare submit command
    # CBIG_CODE_DIR = get_cbig_code_dir()
    CBIG_pbsubmit = os.path.join(CBIG_CODE_DIR, "setup", "CBIG_pbsubmit")
    qsub_command = (
        f"{CBIG_pbsubmit} -cmd \"{cmd_script}\" -walltime {walltime} -mem {memory} -ncpus {ncpus} "
        f"-name {jobname} -joberr {stderr_file} -jobout {stdout_file}")

    # Submit job
    subprocess.run(qsub_command, shell=True, check=True)


def get_number_of_job(username, jobprefix=None):
    """
    Get the total number of submitted jobs for the specified user,
    or the number of jobs matching a specific jobprefix substring.

    Args:
        username (str): Username to query.
        jobprefix (str, optional): Substring to search for in job names.
                                 If None, counts all jobs.

    Returns:
        int: Number of matching jobs (or all jobs if jobprefix is None).
    """

    qstat_command = f"qstat -u {username}"
    result = subprocess.run(qstat_command,
                            shell=True,
                            capture_output=True,
                            text=True)
    lines = result.stdout.strip().split("\n")
    # Exclude the header lines (set as 5 lines for
    # CBIG CSCHPC users)
    job_lines = [line for line in lines[5:] if line.strip()]

    if jobprefix:
        # Count jobs where jobprefix appears in the job name column (
        # set as 4th column (3, by 0 indexing) for CBIG CSCHPC users)
        matching_jobs = [
            line for line in job_lines if jobprefix in line.split()[3]
        ]
        return len(matching_jobs)
    else:
        # If no jobprefix is provided, return the total number of jobs
        return len(job_lines)


def check_and_resubmit_job(check_job_func, submit_job_func, timestamp,
                           try_count, try_max, jobprefix, username,
                           fail_list_flag, check_kwargs, submit_kwargs):
    """Check job status, resubmit jobs if necessary, and retry until success or try_max is reached.

    If fail_list_flag is True, check_job_func is expected to return a
    (success, fail_list) tuple, and fail_list is threaded back into
    check_kwargs/submit_kwargs as 'list_file' for the next attempt (e.g. to
    retry only the failed subset). Otherwise check_job_func is expected to
    return only a success bool.

    Args:
        check_job_func (callable): The function to check if jobs were successful.
        submit_job_func (callable): The function to submit jobs.
        timestamp (str): A string representing the current timestamp
            (accepted by this function but not used in its body).
        try_count (int): The current attempt number for job checking and submission.
        try_max (int): The maximum number of retry attempts.
        jobprefix (str): A prefix for the job name, used in log files and job submission.
        username (str): Username to query in qstat while waiting for jobs to finish.
        fail_list_flag (bool): If True, check_job_func returns a
            (success, fail_list) tuple and fail_list is passed to
            subsequent check/submit calls. If False, check_job_func
            returns only a success bool.
        check_kwargs (dict): Keyword arguments for the check_job_func.
        submit_kwargs (dict): Keyword arguments for the submit_job_func.

    Returns:
        bool: True if all jobs completed successfully (either on the first
            check or after retries). Also returns True (without exhausting
            further retries) if jobprefix contains "NTm" and try_max
            attempts have been exhausted, allowing the caller to proceed
            anyway. Otherwise, the program terminates via sys.exit(1)
            after try_max attempts are exhausted.
    """

    # Initialize the job success status and fail list (if necessary)
    if fail_list_flag:
        all_job_success, fail_list = check_job_func(**check_kwargs)
        submit_kwargs['list_file'] = fail_list
        check_kwargs['list_file'] = fail_list
    else:
        all_job_success = check_job_func(**check_kwargs)

    # If jobs are not successful and max tries haven't been reached, resubmit jobs again
    while not all_job_success and (try_count < (try_max + 1)):
        print(f"Submitting {jobprefix} for try {try_count}/{try_max}...")

        # Submit the job using the provided submit_job_func and arguments
        submit_job_func(**submit_kwargs)

        try_count += 1
        check_kwargs['try_count'] = try_count
        wait_for_job(jobprefix, username)  # Wait for the job to finish

        # Check the job status again
        if fail_list_flag:
            all_job_success, fail_list = check_job_func(**check_kwargs)
            submit_kwargs['list_file'] = fail_list
            check_kwargs['list_file'] = fail_list
        else:
            all_job_success = check_job_func(**check_kwargs)

    # Output final success or failure message
    if all_job_success:
        print(f"{jobprefix} jobs successful.")
        return all_job_success
    elif ("NTm" in jobprefix) and (try_count - 1 == try_max):
        print(f"{jobprefix} jobs failed, attempted {try_count - 1}/{try_max}.")
        print(f"But will still proceed to next selected subID...")
        all_job_success = True
        return all_job_success
    else:
        print(f"{jobprefix} jobs failed, attempted {try_count - 1}/{try_max}.")
        print(f"Terminating script, will not proceed...")
        sys.exit(1)


def run_command_line(command):
    """Run a shell command, terminating the program if it fails.

    Args:
        command (str): The shell command to execute.
    """
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        sys.exit(1)
