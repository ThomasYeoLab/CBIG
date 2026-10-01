'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
"""
This is a library module (no direct sys.argv entry point). It is imported and called by
replication/TMS_Kong2025/script/CBIG_TMSmotion_Kong2025_gt_wrapper.py and
replication/TMS_Kong2025/script/CBIG_TMSmotion_Kong2025_ngt_wrapper.py, e.g.:

Examples:

python $TAN2026_TMSMOTION_DIR/replication/TMS_Kong2025/script/CBIG_TMSmotion_Kong2025_gt_wrapper.py \
    $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json
"""
#############################################################################
# ===========================================================================
# IMPORTS
# ===========================================================================
#############################################################################
import os
import pandas as pd
import sys
import numpy as np
import shutil
import glob
import re

TAN2026_TMSMOTION_DIR = os.environ.get('TAN2026_TMSMOTION_DIR')
sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'util'))

KONG2026_TMS_DIR = os.environ.get('KONG2026_TMS_DIR')
import CBIG_TMSmotion_util as util

CBIG_CODE_DIR = util.get_cbig_code_dir()

#############################################################################
# ===========================================================================
# CONCATENATE WINDOWS TO GENERATE SIMULATION RUNS (ONLY NGT)
# ===========================================================================
#############################################################################


def concat_windows(main_config_file, pl_config_file, run_type, dataset, subID,
                   bld, windowID, window_idx_start, window_idx_end, sim_sesID,
                   sim_bld, preproc_subdir_name, preproc_file_name,
                   out_windowIDs_niis, out_base_dir):
    """Extract one window from a real run's preprocessed NIFTI and append it to a list.

    Uses `fslroi` to extract the frame range [window_idx_start, window_idx_end]
    (inclusive) from the z-scored preprocessed NIFTI of the real run identified by
    `bld`, writes it to a temporary output file, and appends that output file's
    path to `out_windowIDs_niis` so it can later be concatenated (via `fslmerge`)
    with other windows to build a simulated run.

    Args:
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        run_type (str): Run type, e.g. "ngt" or "gt".
        dataset (str): Dataset name.
        subID (str): Subject ID.
        bld (str): Real run (bold) ID that this window is extracted from.
        windowID (str or int): Identifier of this window, used to name the output file.
        window_idx_start (int): Inclusive starting frame index of the window.
        window_idx_end (int): Inclusive ending frame index of the window.
        sim_sesID (str): Simulated session ID this window contributes to.
        sim_bld (str): Simulated bold ID this window contributes to.
        preproc_subdir_name (str): Sub-directory name (relative to the simulated
            run's base dir) that the output window file will live in, expressed
            in terms of `sim_bld`.
        preproc_file_name (str): File name template (expressed in terms of
            `sim_bld`) used to derive both the input real-run file name and the
            output window file name.
        out_windowIDs_niis (list of str): Accumulator list of window NIFTI file
            paths extracted so far; the newly extracted window is appended to it.
        out_base_dir (str): Base output directory for the simulated run.

    Returns:
        list of str: `out_windowIDs_niis` with the newly extracted window's file
        path appended.
    """
    # Read main config
    main_config = util.read_json_file(main_config_file)
    data_dir = main_config['path']['data_dir']

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])

    # Define types
    preproc_type = "after_mc"

    # Modify sub-directory name to appropriate bld/sim_bld
    runID_preproc_subdir_name = preproc_subdir_name.replace(sim_bld, bld)
    sim_runID_preproc_subdir_name = preproc_subdir_name

    # Modify input NIFTI name to appropriate bld/sim_bld
    runID_preproc_file_name = preproc_file_name.replace(
        f'_bld{sim_bld}', f'{subID}_bld{bld}')
    sim_runID_preproc_file_name = preproc_file_name.replace(
        f'_bld{sim_bld}', f'{sim_sesID}_bld{sim_bld}')
    in_runID_nii = os.path.join(data_dir, dataset, 'preprocessed', 'rest',
                                pl_name, f'{preproc_type}_{run_type}', subID,
                                runID_preproc_subdir_name,
                                f'zscore_{runID_preproc_file_name}')

    # check input files exist
    assert os.path.exists(in_runID_nii), f"ERROR: {in_runID_nii} missing!"

    # output dirs
    out_tmp_dir = os.path.join(out_base_dir, sim_runID_preproc_subdir_name,
                               'tmp', sim_runID_preproc_file_name)
    if not os.path.exists(out_tmp_dir):
        os.makedirs(out_tmp_dir)

    # extract z-scored window
    out_windowID_nii = os.path.join(
        out_tmp_dir, f'window_{windowID}_{runID_preproc_file_name}')
    util.run_command_line(
        f'fslroi {in_runID_nii} {out_windowID_nii} {window_idx_start} {window_idx_end - window_idx_start + 1}'
    )

    # append z-scored window
    out_windowIDs_niis.append(out_windowID_nii)

    return out_windowIDs_niis


def extract_append_txt(main_config_file, pl_config_file, run_type, dataset,
                       subID, bld, window_idx_start, window_idx_end,
                       list_subdir_name, list_file_name, sim_sesID, sim_bld,
                       out_windowIDs_lists):
    """Extract one window's rows from a real run's motion-outlier file and append them.

    Loads the real run's motion-outlier (FD or DV) text file, slices out the rows
    corresponding to [window_idx_start, window_idx_end] (inclusive), and appends
    the resulting array to `out_windowIDs_lists` so it can later be stacked with
    other windows to build the motion-outlier file for a simulated run.

    Args:
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        run_type (str): Run type, e.g. "ngt" or "gt".
        dataset (str): Dataset name.
        subID (str): Subject ID.
        bld (str): Real run (bold) ID that this window is extracted from.
        window_idx_start (int): Inclusive starting frame index of the window.
        window_idx_end (int): Inclusive ending frame index of the window.
        list_subdir_name (str): Sub-directory name (relative to the real run's
            preprocessed dir) containing the motion-outlier file.
        list_file_name (str): File name template (expressed in terms of
            `sim_bld`) used to derive the real-run motion-outlier file name.
        sim_sesID (str): Simulated session ID this window contributes to.
        sim_bld (str): Simulated bold ID this window contributes to.
        out_windowIDs_lists (list of np.ndarray): Accumulator list of
            window-level motion-outlier arrays extracted so far; the newly
            extracted window's array is appended to it.

    Returns:
        list of np.ndarray: `out_windowIDs_lists` with the newly extracted
        window's array appended.
    """
    # Read main config
    main_config = util.read_json_file(main_config_file)
    data_dir = main_config['path']['data_dir']

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])

    # Define types
    preproc_type = "after_mc"

    # input files
    runID_list_file_name = list_file_name.replace(f'_bld{sim_bld}',
                                                  f'{subID}_bld{bld}')
    in_runID_outliers = os.path.join(data_dir, dataset, 'preprocessed', 'rest',
                                     pl_name, f'{preproc_type}_{run_type}',
                                     subID, list_subdir_name,
                                     runID_list_file_name)

    # check input files exist
    assert os.path.exists(
        in_runID_outliers), f"ERROR: {in_runID_outliers} missing!"
    # skip reference frame row
    arr_runID_outliers = np.loadtxt(in_runID_outliers, dtype=float)
    arr_runID_outliers = np.reshape(arr_runID_outliers, (-1, 1))
    arr_windowID_outliers = arr_runID_outliers[
        window_idx_start:window_idx_end + 1]
    out_windowIDs_lists.append(arr_windowID_outliers)

    return out_windowIDs_lists


def convert_list2array(main_config_file, pl_config_file, dataset, subID,
                       sim_run_bin, sim_sesID, sim_bld, list_subdir_name,
                       list_file_name, expected_frame_count, out_base_dir,
                       out_windowIDs_lists):
    """Stack accumulated window-level motion-outlier arrays and save as a text file.

    Vertically stacks the list of per-window motion-outlier arrays into a single
    array for the simulated run, asserts that the resulting frame count matches
    the expected total, and saves it as a text file under `out_base_dir`.

    Args:
        main_config_file (str): Path to the main pipeline config JSON file
            (read but unused beyond validating the file exists).
        pl_config_file (str): Path to the pipeline-specific config JSON file
            (read but unused beyond validating the file exists).
        dataset (str): Dataset name (used only in the assertion error message).
        subID (str): Subject ID (used only in the assertion error message).
        sim_run_bin (str): Simulation run bin identifier (used only in the
            assertion error message).
        sim_sesID (str): Simulated session ID used to derive the output file name.
        sim_bld (str): Simulated bold ID used to derive the output file name.
        list_subdir_name (str): Sub-directory (relative to `out_base_dir`) that
            the output file is saved into.
        list_file_name (str): File name template (expressed in terms of
            `sim_bld`) used to derive the output file name.
        expected_frame_count (int): Expected number of frames (rows) in the
            stacked array; an assertion fails if the actual count differs.
        out_base_dir (str): Base output directory for the simulated run.
        out_windowIDs_lists (list of np.ndarray): Per-window motion-outlier
            arrays to stack together.
    """
    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    util.read_json_file(pl_config_file)

    # input files
    sim_runID_list_file_name = list_file_name.replace(
        f'_bld{sim_bld}', f'{sim_sesID}_bld{sim_bld}')

    # output dirs
    out_dir = os.path.join(out_base_dir, list_subdir_name)
    os.makedirs(out_dir, exist_ok=True)

    # convert to np array
    out_windowIDs_lists = np.vstack(out_windowIDs_lists)

    # check expected number of frames
    assert len(out_windowIDs_lists) == expected_frame_count, (
        f"ERROR: Expected {expected_frame_count}, "
        f"but got {len(out_windowIDs_lists)} frames "
        f"for {dataset}, {subID}, {sim_run_bin}, {sim_sesID}, {sim_bld}.")

    # save txt file
    np.savetxt(os.path.join(out_dir, sim_runID_list_file_name),
               out_windowIDs_lists,
               fmt='%.6f')


def concat_wrapper(main_config_file, pl_config_file, BOLD_stem,
                   BOLD_stem_pre_mc, dataset, subID, sim_run_bin, sim_sesID,
                   sim_bld, sesID, runID, window_list, run_type,
                   reduce_mem_flag):
    """Build one simulated run by concatenating its constituent windows.

    Top-level wrapper (intended to be run inside a submitted PBS job via
    `submit_concat_job`) that assembles a single simulated run (`sim_sesID`,
    `sim_bld`) out of the windows listed in `window_list`. For each window row
    it extracts the corresponding volume/surface NIFTI segments (via
    `concat_windows`) and motion-outlier rows (via `extract_append_txt`),
    accumulating them across windows. It then merges the NIFTI segments with
    `fslmerge`, removes temporary per-window files/dirs, and writes out the
    concatenated motion-outlier files (via `convert_list2array`).

    Args:
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        BOLD_stem (str): BOLD preprocessing stem used to build preprocessed
            file names.
        BOLD_stem_pre_mc (str): BOLD preprocessing stem (pre motion-correction)
            used to build motion-outlier file names.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        sim_run_bin (str): Simulation run bin identifier.
        sim_sesID (str): Simulated session ID being built.
        sim_bld (str): Simulated bold ID being built.
        sesID (str): Real session ID (passed through, not directly used here).
        runID (str): Real run ID (passed through, not directly used here).
        window_list (str): Path to a CSV file listing the windows (with columns
            such as `windowID`, `window_idx_start`, `window_idx_end`,
            `sim_sesID`, `bld`, `sim_bld`, `TR`) that make up this simulated run.
        run_type (str): Run type, e.g. "ngt".
        reduce_mem_flag (str): "1" to use the reduced-memory (loose-masked,
            non-T1-space) preprocessed volume file, or "0" to use the full
            T1-space volume file.
    """
    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    func2T1_sm = pl_config['func2T1']['func2T1_sm']
    window_time_sec = pl_config['select_subID']['window']['window_time_sec']

    # Init tolerance for assert statements
    tol_window_duration = 5

    # Assign input dirs
    sim_ses_preproc_real_dir = util.get_sim_ses_preproc_dir_real_path(
        main_config_file, pl_config_file, dataset, subID, sim_run_bin,
        sim_sesID)

    # initialize NIFTI lists
    dict_preproc = {}
    if reduce_mem_flag == "1":
        dict_key = f'_bld{sim_bld}{BOLD_stem}_sm{func2T1_sm}_loosemask.nii.gz'
        dict_preproc[dict_key] = {os.path.join('bold', sim_bld): []}
    elif reduce_mem_flag == "0":
        dict_key = f'_bld{sim_bld}{BOLD_stem}_T1_sm{func2T1_sm}_loosemask.nii.gz'
        dict_preproc[dict_key] = {os.path.join('bold', sim_bld): []}
    else:
        raise ValueError(
            f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
        )
    dict_preproc[f'lh._bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz'] = {
        os.path.join('surf'): []
    }
    dict_preproc[f'rh._bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz'] = {
        os.path.join('surf'): []
    }

    # initialize txt lists
    dict_lists = {}
    dict_lists[f'_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{FD_stem}'] = {
        os.path.join('bold', 'mc'): []
    }
    dict_lists[f'_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{DV_stem}'] = {
        os.path.join('bold', 'mc'): []
    }

    # Init counter
    expected_frame_count = 0
    # Loop through each windowID
    window_df = pd.read_csv(window_list, dtype={'bld': str, 'sim_bld': str})
    for idx_windowID, row_windowID in window_df.iterrows():
        # identifier variables
        windowID = row_windowID['windowID']
        window_idx_start = row_windowID['window_idx_start']  # inclusive
        window_idx_end = row_windowID['window_idx_end']  # inclusive
        window_sim_sesID = row_windowID['sim_sesID']
        window_bld = row_windowID['bld']
        window_sim_bld = row_windowID['sim_bld']
        window_TR = row_windowID['TR']

        # Check identities
        assert window_sim_sesID == sim_sesID, (
            f"ERROR: {dataset} {subID} window_sim_sesID "
            f"({window_sim_sesID}) != sim_sesID ({sim_sesID}).")
        assert window_sim_bld == sim_bld, (
            f"ERROR: {dataset} {subID} window_sim_bld "
            f"({window_sim_bld}) != sim_bld ({sim_bld}).")

        # Check windows of expected duration
        assert abs((window_TR * (window_idx_end - window_idx_start + 1)) - window_time_sec) < tol_window_duration, \
            (f'ERROR: window duration != window_time_sec for {dataset}, {subID}, {sim_run_bin}, '
             f'{window_sim_sesID}, {window_sim_bld}, {windowID}')

        # generate nii file
        # concat z-scored windows
        for preproc_file_name in dict_preproc.keys():
            preproc_subdir_name = next(iter(dict_preproc[preproc_file_name]))
            dict_preproc[preproc_file_name][
                preproc_subdir_name] = concat_windows(
                    main_config_file, pl_config_file, run_type, dataset, subID,
                    window_bld, windowID, window_idx_start, window_idx_end,
                    window_sim_sesID, window_sim_bld, preproc_subdir_name,
                    preproc_file_name,
                    dict_preproc[preproc_file_name][preproc_subdir_name],
                    sim_ses_preproc_real_dir)

        # generate motion_outlier files
        # extract window-relevant data & append
        for list_file_name in dict_lists.keys():
            list_subdir_name = next(iter(dict_lists[list_file_name]))
            out_windowIDs_lists = dict_lists[list_file_name][list_subdir_name]
            out_windowIDs_lists = extract_append_txt(
                main_config_file, pl_config_file, run_type, dataset, subID,
                window_bld, window_idx_start, window_idx_end, list_subdir_name,
                list_file_name, window_sim_sesID, window_sim_bld,
                out_windowIDs_lists)
        # Accumulate total frame count for simulation run
        # for later assert statement
        expected_frame_count += (window_idx_end - window_idx_start + 1)

    # generate nii file
    for preproc_file_name in dict_preproc.keys():
        sim_runID_preproc_file_name = preproc_file_name.replace(
            f'_bld{sim_bld}', f'{sim_sesID}_bld{sim_bld}')
        # Modify sub-directory name to appropriate bld/sim_bld
        preproc_subdir_name = next(iter(dict_preproc[preproc_file_name]))
        # concatenate windows
        merged_output = os.path.join(sim_ses_preproc_real_dir,
                                     preproc_subdir_name,
                                     sim_runID_preproc_file_name)
        util.run_command_line(
            f'fslmerge -t {merged_output} {" ".join(dict_preproc[preproc_file_name][preproc_subdir_name])} 2>/dev/null'
        )

        # delete tmp subdir
        shutil.rmtree(
            os.path.join(sim_ses_preproc_real_dir, preproc_subdir_name, 'tmp',
                         sim_runID_preproc_file_name))

    # delete tmp dir
    patterns = [os.path.join(sim_ses_preproc_real_dir, 'bold', sim_bld, 'tmp')]
    for pattern in patterns:
        for tmp_dir in glob.glob(pattern):
            if os.path.isdir(tmp_dir):
                shutil.rmtree(tmp_dir)

    # generate motion_outlier files
    # convert motion_outliers lists to np arrays
    for list_file_name in dict_lists.keys():
        list_subdir_name = next(iter(dict_lists[list_file_name]))
        out_windowIDs_lists = dict_lists[list_file_name][list_subdir_name]
        convert_list2array(main_config_file, pl_config_file, dataset, subID,
                           sim_run_bin, sim_sesID, sim_bld, list_subdir_name,
                           list_file_name, expected_frame_count,
                           sim_ses_preproc_real_dir, out_windowIDs_lists)


def submit_concat_job(timestamp, jobprefix, main_config_file, pl_config_file,
                      dataset, subID, sim_run_bin, sim_sesID, sim_bld, sesID,
                      runID, window_list, run_type, reduce_mem_flag):
    """Submit a PBS job that builds one simulated run via `concat_wrapper`.

    Deletes any stale outputs from a previous attempt (the simulated bold dir,
    outlier files, and surface files), then submits a job that re-invokes this
    module's `concat_wrapper` function with the given parameters. Memory and
    walltime are chosen based on `dataset`/`subID` and `reduce_mem_flag`.

    Args:
        timestamp (str): Timestamp identifying this pipeline run, used to
            namespace job output directories.
        jobprefix (str): Prefix used to build the PBS job name.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        sim_run_bin (str): Simulation run bin identifier.
        sim_sesID (str): Simulated session ID being built.
        sim_bld (str): Simulated bold ID being built.
        sesID (str): Real session ID (passed through to `concat_wrapper`).
        runID (str): Real run ID (passed through to `concat_wrapper`).
        window_list (str): Path to the CSV file listing the windows that make
            up this simulated run.
        run_type (str): Run type; must be "ngt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.
    """

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    if run_type == "ngt":
        rep_TMS_Kong2025_concat_dir = pl_config['path'][
            'rep_TMS_Kong2025_ngt_concat_dir']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'ngt'.")

    # Define types
    preproc_type = "after_mc"

    # Determine file paths and directories
    resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                       pl_config_file,
                                                       timestamp,
                                                       preproc_type,
                                                       dataset,
                                                       subID,
                                                       run_type=run_type)
    std_dir = os.path.join(rep_TMS_Kong2025_concat_dir, timestamp, "job")
    os.makedirs(std_dir, exist_ok=True)

    # Initialize input vars
    BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                  resume_config_file)
    BOLD_stem_pre_mc = util.generate_pre_mc_stem_from_precursor_stem(BOLD_stem)

    # Clean up old output(s)
    # Delete `sim_bld` bold dir
    sim_ses_preproc_real_dir = util.get_sim_ses_preproc_dir_real_path(
        main_config_file, pl_config_file, dataset, subID, sim_run_bin,
        sim_sesID)
    sim_bld_preproc_real_dir = os.path.join(sim_ses_preproc_real_dir, "bold",
                                            sim_bld)
    util.restart_dir(sim_bld_preproc_real_dir)
    # Delete outliers files
    sim_bld_mc_preproc_real_dir = os.path.join(sim_ses_preproc_real_dir,
                                               "bold", "mc")
    FD_metric_outliers_file = os.path.join(
        sim_bld_mc_preproc_real_dir,
        f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{FD_stem}")
    DV_metric_outliers_file = os.path.join(
        sim_bld_mc_preproc_real_dir,
        f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{DV_stem}")
    util.delete_file(FD_metric_outliers_file)
    util.delete_file(DV_metric_outliers_file)
    # Delete surface files
    sim_bld_surf_preproc_real_dir = os.path.join(sim_ses_preproc_real_dir,
                                                 "surf")
    lh_file = os.path.join(
        sim_bld_surf_preproc_real_dir,
        f'lh.{sim_sesID}_bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz')
    rh_file = os.path.join(
        sim_bld_surf_preproc_real_dir,
        f'rh.{sim_sesID}_bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz')
    util.delete_file(lh_file)
    util.delete_file(rh_file)

    # PBS job parameters
    jobname = f"{jobprefix}_{dataset}_{subID}_{sim_run_bin}_{sim_sesID}_{sim_bld}"
    stderr_file = os.path.join(
        std_dir, f"{dataset}_{subID}_{sim_run_bin}_{sim_sesID}_{sim_bld}.err")
    stdout_file = os.path.join(
        std_dir, f"{dataset}_{subID}_{sim_run_bin}_{sim_sesID}_{sim_bld}.out")

    # Directly construct the command script for job submission
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; "
        f"sys.path.append(\\\"{TAN2026_TMSMOTION_DIR}/TMS/Kong2025\\\"); "
        f"import CBIG_TMSmotion_Kong2025 as TMS_Kong2025; "
        f"TMS_Kong2025.concat_wrapper("
        f"\\\"{main_config_file}\\\", \\\"{pl_config_file}\\\", \\\"{BOLD_stem}\\\", "
        f"\\\"{BOLD_stem_pre_mc}\\\", \\\"{dataset}\\\", \\\"{subID}\\\", \\\"{sim_run_bin}\\\", "
        f"\\\"{sim_sesID}\\\", \\\"{sim_bld}\\\", \\\"{sesID}\\\", \\\"{runID}\\\", "
        f"\\\"{window_list}\\\", \\\"{run_type}\\\", \\\"{reduce_mem_flag}\\\");' "
        f"source deactivate")

    # Define memory & walltime
    if reduce_mem_flag == "1":
        walltime = "1:00:00"
        very_high_memory = "8G"
        high_memory = "8G"
        intermediate_memory = "8G"
        normal_memory = "8G"
        if dataset == "Newbold" and (subID == "sub-cast2"
                                     or subID == "sub-cast3"):
            memory = very_high_memory
        elif dataset == "MSC" or (dataset == "Newbold"
                                  and subID == "sub-cast1"):
            memory = high_memory
        elif dataset == "Myconnectome":
            memory = intermediate_memory
        else:
            memory = normal_memory
    elif reduce_mem_flag == "0":
        walltime = "3:00:00"
        very_high_memory = "240G"
        high_memory = "130G"
        intermediate_memory = "100G"
        normal_memory = "60G"
        if dataset == "Newbold" and (subID == "sub-cast2"
                                     or subID == "sub-cast3"):
            memory = very_high_memory
        elif dataset == "MSC" or (dataset == "Newbold"
                                  and subID == "sub-cast1"):
            memory = high_memory
        elif dataset == "Myconnectome":
            memory = intermediate_memory
        else:
            memory = normal_memory
    else:
        raise ValueError(
            f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
        )

    # Submit job
    util.submit_job(walltime=walltime,
                    memory=memory,
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def submit_all_concat_job(timestamp, jobprefix, main_config_file,
                          pl_config_file, run_type, list_file,
                          reduce_mem_flag):
    """Submit concat jobs for every simulation run listed in `list_file`.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        jobprefix (str): Prefix used to build PBS job names.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        run_type (str): Run type; must be "ngt".
        list_file (str): Path to a CSV file listing simulation runs (one row
            per `dataset`/`subID`/`sim_run_bin`/`sim_sesID`/`sim_bld`) to
            submit concat jobs for.
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']

    # Read list of simulation runs across all subIDs
    sim_run_df = pd.read_csv(list_file, dtype={'bld': str, 'sim_bld': str})
    # Loop through each simulation run
    for idx, row in sim_run_df.iterrows():
        dataset, subID, sim_run_bin, sim_sesID, sim_bld, sesID, runID = row[[
            'dataset', 'subID', 'sim_run_bin', 'sim_sesID', 'sim_bld', 'sesID',
            'runID'
        ]]

        # For that 1 sim_bld, read list of windows belonging to that sim_bld
        window_list = os.path.join(rep_sel_subID_sim_dir, 'result', 'summary',
                                   'window_list_1', dataset, subID,
                                   sim_run_bin, sim_sesID, sim_bld,
                                   'window_list_1.txt')

        # Submit job for 1 simulation run
        submit_concat_job(timestamp, jobprefix, main_config_file,
                          pl_config_file, dataset, subID, sim_run_bin,
                          sim_sesID, sim_bld, sesID, runID, window_list,
                          run_type, reduce_mem_flag)


def check_concat_job(timestamp, try_count, main_config_file, pl_config_file,
                     data_dir, dataset, subID, sim_run_bin, sim_sesID, sim_bld,
                     sesID, runID, all_job_success, run_type, reduce_mem_flag):
    """Check whether one concat job's output files were created successfully.

    Verifies that the surface, bold, and motion-outlier files expected from a
    `submit_concat_job` run exist. Any missing file is logged to a per-attempt
    failure message file, and the corresponding run identifiers are logged to a
    per-attempt failure list file.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number, used to namespace the failure log files.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        data_dir (str): Data directory (read but unused beyond being passed
            through).
        dataset (str): Dataset name.
        subID (str): Subject ID.
        sim_run_bin (str): Simulation run bin identifier.
        sim_sesID (str): Simulated session ID.
        sim_bld (str): Simulated bold ID.
        sesID (str): Real session ID (used only in failure log messages).
        runID (str): Real run ID (used only in failure log messages).
        all_job_success (bool): Running success flag accumulated across calls;
            set to False if this run's files are missing.
        run_type (str): Run type; must be "ngt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.

    Returns:
        tuple: (bool, str) the updated `all_job_success` flag, and the path to
        the failure list file for this `try_count`.
    """

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    func2T1_sm = pl_config['func2T1']['func2T1_sm']
    if run_type == "ngt":
        rep_TMS_Kong2025_concat_dir = pl_config['path'][
            'rep_TMS_Kong2025_ngt_concat_dir']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'ngt'.")

    # Define types
    preproc_type = "after_mc"

    # Define input dirs
    sim_ses_preproc_real_dir = util.get_sim_ses_preproc_dir_real_path(
        main_config_file, pl_config_file, dataset, subID, sim_run_bin,
        sim_sesID)

    # Define files to check
    resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                       pl_config_file,
                                                       timestamp,
                                                       preproc_type,
                                                       dataset,
                                                       subID,
                                                       run_type=run_type)
    BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                  resume_config_file)
    BOLD_stem_pre_mc = util.generate_pre_mc_stem_from_precursor_stem(BOLD_stem)

    filename2check = [
        f"surf/lh.{sim_sesID}_bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz",
        f"surf/rh.{sim_sesID}_bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz",
        f"bold/mc/{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{FD_stem}",
        f"bold/mc/{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{DV_stem}"
    ]
    if reduce_mem_flag == "1":
        filename2check.append(
            f"bold/{sim_bld}/{sim_sesID}_bld{sim_bld}{BOLD_stem}_sm{func2T1_sm}_loosemask.nii.gz"
        )
    elif reduce_mem_flag == "0":
        filename2check.append(
            f"bold/{sim_bld}/{sim_sesID}_bld{sim_bld}{BOLD_stem}_T1_sm{func2T1_sm}_loosemask.nii.gz"
        )
    else:
        raise ValueError(
            f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
        )

    # Define output files
    timestamp_dir = os.path.join(rep_TMS_Kong2025_concat_dir, timestamp,
                                 f"{dataset}_{subID}")
    os.makedirs(timestamp_dir, exist_ok=True)
    fail_list = os.path.join(timestamp_dir, f"fail_list_try_{try_count}.txt")
    fail_msg = os.path.join(timestamp_dir, f"fail_msg_try_{try_count}.txt")

    # Check job output files for success
    id_success = True
    for filename in filename2check:
        exist = util.check_job_output(
            os.path.join(sim_ses_preproc_real_dir, filename))
        if not exist:
            id_success = False
            # Handle failure message list
            if not os.path.exists(fail_msg):
                with open(fail_msg, 'w') as f:
                    f.write(
                        "dataset,subID,sim_run_bin,sim_sesID,sim_bld,sesID,runID    error_message\n"
                    )
            with open(fail_msg, 'a') as f:
                f.write(
                    f"{dataset},{subID},{sim_run_bin},{sim_sesID},{sim_bld},{sesID},{runID}    {filename}\n"
                )
    if not id_success:
        all_job_success = False
        # Handle failed ID list
        if not os.path.exists(fail_list):
            with open(fail_list, 'w') as f:
                f.write(
                    "dataset,subID,sim_run_bin,sim_sesID,sim_bld,sesID,runID\n"
                )
        with open(fail_list, 'a') as f:
            f.write(
                f"{dataset},{subID},{sim_run_bin},{sim_sesID},{sim_bld},{sesID},{runID}\n"
            )

    return all_job_success, fail_list


def check_all_concat_job(timestamp, try_count, main_config_file,
                         pl_config_file, data_dir, list_file, run_type,
                         reduce_mem_flag):
    """Check every concat job listed in `list_file` for successful completion.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number, used to namespace the failure log files.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        data_dir (str): Data directory, passed through to `check_concat_job`.
        list_file (str): Path to a CSV file listing simulation runs to check.
        run_type (str): Run type; must be "ngt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.

    Returns:
        tuple: (bool, str) whether all checked jobs succeeded, and the path to
        the failure list file for this `try_count` (None if `list_file` is empty).
    """

    # Read list of simulation runs across all subIDs
    sim_run_df = pd.read_csv(list_file, dtype={'bld': str, 'sim_bld': str})

    # Init check job vars
    all_job_success = True
    fail_list = None

    # Loop through each simulation run
    for idx, row in sim_run_df.iterrows():
        dataset, subID, sim_run_bin, sim_sesID, sim_bld, sesID, runID = row[[
            'dataset', 'subID', 'sim_run_bin', 'sim_sesID', 'sim_bld', 'sesID',
            'runID'
        ]]

        # Check job
        all_job_success, fail_list = check_concat_job(
            timestamp, try_count, main_config_file, pl_config_file, data_dir,
            dataset, subID, sim_run_bin, sim_sesID, sim_bld, sesID, runID,
            all_job_success, run_type, reduce_mem_flag)

    return all_job_success, fail_list


#############################################################################
# ===========================================================================
# GENERATE MOTION OUTLIERS FILES (ONLY NGT)
# ===========================================================================
#############################################################################
def submit_outlier_job(timestamp, jobprefix, main_config_file, pl_config_file,
                       dataset, subID, sim_run_bin, sim_sesID, sim_sesID_df,
                       run_type):
    """Submit one PBS job that computes motion outliers for all thresholds of a sim session.

    Builds a single MATLAB job (calling `CBIG_TMS_motion_outliers` once
    per row of `sim_sesID_df`) that computes the FD/DV motion-outlier file for
    every varied FD/DV threshold combination belonging to this simulated
    session (`sim_sesID`), so that all threshold combinations for one session
    share a single submitted job. Deletes any stale outlier file before
    appending each threshold's MATLAB command.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        jobprefix (str): Prefix used to build the PBS job name.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        sim_run_bin (str): Simulation run bin identifier.
        sim_sesID (str): Simulated session ID.
        sim_sesID_df (pd.DataFrame): Rows (one per `sim_bld`/threshold
            combination) belonging to this simulated session, with columns
            including `sim_bld`, `vary_<FD_metric>_th`, and `vary_<DV_metric>_th`.
        run_type (str): Run type; must be "ngt".
    """

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rm_seg = pl_config['preproc']['rm_seg']
    oneBtwoA_flag = pl_config['TMS']['oneBtwoA_flag']
    min5_flag = pl_config['TMS']['min5_flag']
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    if run_type == "ngt":
        rep_TMS_Kong2025_outlier_dir = pl_config['path'][
            'rep_TMS_Kong2025_ngt_outlier_dir']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'ngt'.")

    # Define types
    preproc_type = "after_mc"

    # Initialize MATLAB command
    matlab_cmds = ["clear; clc; close all;"]
    matlab_cmds.append(f"addpath('{KONG2026_TMS_DIR}/targeting');")

    # Loop through each varied FD & DV threshold set
    # to keep appending MATLAB command to submit
    # multiple threshold sets in a single job
    for idx, row in sim_sesID_df.iterrows():
        sim_bld = row['sim_bld']
        vary_FD_th = row[f"vary_{FD_metric}_th"]
        vary_DV_th = row[f"vary_{DV_metric}_th"]

        # Determine file paths and directories
        resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                           pl_config_file,
                                                           timestamp,
                                                           preproc_type,
                                                           dataset,
                                                           subID,
                                                           run_type=run_type)
        std_dir = os.path.join(rep_TMS_Kong2025_outlier_dir, timestamp, "job",
                               dataset, subID, sim_run_bin, sim_sesID, sim_bld)
        os.makedirs(std_dir, exist_ok=True)
        sim_run_preproc_dir = util.get_sim_ses_preproc_dir_real_path(
            main_config_file, pl_config_file, dataset, subID, sim_run_bin,
            sim_sesID)

        # Initialize input vars
        BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                      resume_config_file)
        BOLD_stem_pre_mc = util.generate_pre_mc_stem_from_precursor_stem(
            BOLD_stem)
        DV_outlier_file = os.path.join(
            sim_run_preproc_dir, "bold", "mc",
            f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{DV_stem}"
        )
        FD_outlier_file = os.path.join(
            sim_run_preproc_dir, "bold", "mc",
            f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{FD_stem}"
        )
        qc_dir = os.path.join(sim_run_preproc_dir, "qc")
        output_dir_name = os.path.join(qc_dir, f"{sim_sesID}_bld{sim_bld}")
        os.makedirs(qc_dir, exist_ok=True)

        # Clean up old output(s)
        outlier_file = os.path.join(
            qc_dir,
            f"{sim_sesID}_bld{sim_bld}{FD_stem}{vary_FD_th}{DV_stem}{vary_DV_th}_motion_outliers.txt"
        )
        util.delete_file(outlier_file)

        # Append MATLAB command for current varied threshold
        matlab_cmds.append("clear; clc; close all;")
        matlab_cmds.append(
            f"CBIG_TMS_motion_outliers(\'{FD_stem}\', \'{DV_outlier_file}\', \'{FD_outlier_file}\', "
            f"\'{vary_FD_th}\', \'{vary_DV_th}\', \'{rm_seg}\', \'{output_dir_name}\', "
            f"\'{oneBtwoA_flag}\', \'{min5_flag}\');")

    # PBS job parameters
    jobname = f"{jobprefix}_{dataset}_{subID}_{sim_run_bin}_{sim_sesID}"
    stderr_file = os.path.join(
        std_dir, f"{dataset}_{subID}_{sim_run_bin}_{sim_sesID}.err")
    stdout_file = os.path.join(
        std_dir, f"{dataset}_{subID}_{sim_run_bin}_{sim_sesID}.out")

    # Directly construct the command script for job submission
    matlab_cmds.append("exit;")
    full_matlab_cmd = " ".join(matlab_cmds)
    cmd_script = f"matlab -nosplash -nodisplay -nodesktop -r \\\"{full_matlab_cmd}\\\""

    # Submit job
    util.submit_job(walltime="1:00:00",
                    memory="2G",
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def submit_all_outlier_job(timestamp, jobprefix, main_config_file,
                           pl_config_file, run_type, list_file):
    """Submit outlier jobs for every simulated session listed in `list_file`.

    Groups `list_file` by simulated session (`dataset`, `subID`, `sim_run_bin`,
    `sim_sesID`) and submits one job per group via `submit_outlier_job`.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        jobprefix (str): Prefix used to build PBS job names.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        run_type (str): Run type; must be "ngt".
        list_file (str): Path to a CSV file listing simulation runs and their
            varied FD/DV thresholds.
    """
    # Read pipeline config
    util.read_json_file(pl_config_file)

    # Read list of simulation runs across all thresholds across all subIDs
    sim_run_df = pd.read_csv(list_file, dtype={'bld': str, 'sim_bld': str})

    # Group by simulation session ID-level
    group_cols = ['dataset', 'subID', 'sim_run_bin', 'sim_sesID']
    grouped = sim_run_df.groupby(group_cols)

    # Loop through each simulation session ID
    for group_keys, sim_sesID_df in grouped:
        dataset, subID, sim_run_bin, sim_sesID = group_keys

        # Submit 1 job for each simulation session ID
        submit_outlier_job(timestamp, jobprefix, main_config_file,
                           pl_config_file, dataset, subID, sim_run_bin,
                           sim_sesID, sim_sesID_df, run_type)


def check_outlier_job(timestamp, try_count, main_config_file, pl_config_file,
                      data_dir, dataset, subID, sim_run_bin, sim_sesID,
                      sim_bld, vary_FD_th, vary_DV_th, sesID, runID,
                      all_job_success, run_type):
    """Check whether one threshold combination's outlier file was created successfully.

    Verifies that the motion-outlier text file for one `sim_bld` and one
    FD/DV threshold combination exists. Any missing file is logged to a
    per-attempt failure message file, and the corresponding identifiers are
    logged to a per-attempt failure list file.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number, used to namespace the failure log files.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        data_dir (str): Data directory (read but unused beyond being passed
            through).
        dataset (str): Dataset name.
        subID (str): Subject ID.
        sim_run_bin (str): Simulation run bin identifier.
        sim_sesID (str): Simulated session ID.
        sim_bld (str): Simulated bold ID.
        vary_FD_th (float or str): Varied FD threshold for this combination.
        vary_DV_th (float or str): Varied DV threshold for this combination.
        sesID (str): Real session ID (used only in failure log messages).
        runID (str): Real run ID (used only in failure log messages).
        all_job_success (bool): Running success flag accumulated across calls;
            set to False if this combination's file is missing.
        run_type (str): Run type; must be "ngt".

    Returns:
        tuple: (bool, str) the updated `all_job_success` flag, and the path to
        the failure list file for this `try_count`.
    """

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    if run_type == "ngt":
        rep_TMS_Kong2025_outlier_dir = pl_config['path'][
            'rep_TMS_Kong2025_ngt_outlier_dir']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'ngt'.")

    # Define files to check
    filename2check = [
        f"{sim_sesID}_bld{sim_bld}{FD_stem}{vary_FD_th}{DV_stem}{vary_DV_th}_motion_outliers.txt"
    ]
    sim_run_preproc_dir = util.get_sim_ses_preproc_dir_real_path(
        main_config_file, pl_config_file, dataset, subID, sim_run_bin,
        sim_sesID)
    qc_dir = os.path.join(sim_run_preproc_dir, "qc")

    # Define output files
    timestamp_dir = os.path.join(rep_TMS_Kong2025_outlier_dir, timestamp,
                                 f"{dataset}_{subID}")
    os.makedirs(timestamp_dir, exist_ok=True)
    fail_list = os.path.join(timestamp_dir, f"fail_list_try_{try_count}.txt")
    fail_msg = os.path.join(timestamp_dir, f"fail_msg_try_{try_count}.txt")

    # Check job output files for success
    id_success = True
    for filename in filename2check:
        exist = util.check_job_output(os.path.join(qc_dir, filename))
        if not exist:
            id_success = False
            # Handle failure message list
            if not os.path.exists(fail_msg):
                with open(fail_msg, 'w') as f:
                    f.write(
                        f"dataset,subID,sim_run_bin,sim_sesID,sim_bld,vary_{FD_metric}_th,"
                        f"vary_{DV_metric}_th,sesID,runID    error_message\n")
            with open(fail_msg, 'a') as f:
                f.write(
                    f"{dataset},{subID},{sim_run_bin},{sim_sesID},{sim_bld},{vary_FD_th},"
                    f"{vary_DV_th},{sesID},{runID}    {filename}\n")
    if not id_success:
        all_job_success = False
        # Handle failed ID list
        if not os.path.exists(fail_list):
            with open(fail_list, 'w') as f:
                f.write(
                    f"dataset,subID,sim_run_bin,sim_sesID,sim_bld,vary_{FD_metric}_th,vary_{DV_metric}_th,sesID,runID\n"
                )
        with open(fail_list, 'a') as f:
            f.write(
                f"{dataset},{subID},{sim_run_bin},{sim_sesID},{sim_bld},{vary_FD_th},{vary_DV_th},{sesID},{runID}\n"
            )

    return all_job_success, fail_list


def check_all_outlier_job(timestamp, try_count, main_config_file,
                          pl_config_file, data_dir, list_file, run_type):
    """Check every outlier job (per simulation run and threshold) in `list_file`.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number, used to namespace the failure log files.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        data_dir (str): Data directory, passed through to `check_outlier_job`.
        list_file (str): Path to a CSV file listing simulation runs and their
            varied FD/DV thresholds.
        run_type (str): Run type; must be "ngt".

    Returns:
        tuple: (bool, str) whether all checked jobs succeeded, and the path to
        the failure list file for this `try_count` (None if `list_file` is empty).
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']

    # Read list of simulation runs across all subIDs
    sim_run_df = pd.read_csv(list_file, dtype={'bld': str, 'sim_bld': str})

    # Init check job vars
    all_job_success = True
    fail_list = None

    # Loop through each simulation run for each threshold
    for idx, row in sim_run_df.iterrows():
        dataset, subID, sim_run_bin, sim_sesID, sim_bld, vary_FD_th, vary_DV_th, sesID, runID = row[
            [
                'dataset', 'subID', 'sim_run_bin', 'sim_sesID', 'sim_bld',
                f"vary_{FD_metric}_th", f"vary_{DV_metric}_th", 'sesID',
                'runID'
            ]]

        # Check job
        all_job_success, fail_list = check_outlier_job(
            timestamp, try_count, main_config_file, pl_config_file, data_dir,
            dataset, subID, sim_run_bin, sim_sesID, sim_bld, vary_FD_th,
            vary_DV_th, sesID, runID, all_job_success, run_type)

    return all_job_success, fail_list


#############################################################################
# ===========================================================================
# GENERATE SIMULATION DIR STORING SYMBOLIC LINKS (ONLY NGT)
# ===========================================================================
#############################################################################
def generate_symlink(file_subdir_dict,
                     base_src_dir,
                     base_tgt_dir,
                     bld=None,
                     new_bld=None,
                     subID=None,
                     sim_sesID=None):
    """Create symbolic links from a source to a target "CBIG preproc structured" dir.

    For each file in `file_subdir_dict`, creates a symbolic link under
    `base_tgt_dir` pointing to the corresponding real file under
    `base_src_dir`. Optionally renames the sub-directory/file name in the
    target by substituting `bld` -> `new_bld` and `subID` -> `sim_sesID`
    (while guarding against accidentally corrupting a `sub-XXX`-style ID
    that happens to contain the `bld` substring). After each link is
    created, asserts that it exists and correctly resolves, which avoids
    the need for a separate `check_XXX` function for this step.

    Args:
        file_subdir_dict (dict): Mapping of {file_name: sub_dir_name}
            identifying, for each file to link, the sub-directory (relative
            to `base_src_dir`/`base_tgt_dir`) it lives in.
        base_src_dir (str): Base directory containing the real source files.
        base_tgt_dir (str): Base directory in which the symlinks are created.
        bld (str, optional): Source bold ID to replace in the target
            sub-directory/file name. Defaults to None (no bld substitution).
        new_bld (str, optional): Replacement bold ID substituted for `bld` in
            the target sub-directory/file name. Defaults to None.
        subID (str, optional): Source subject ID to replace with `sim_sesID`
            in the target file name. Defaults to None (no subID substitution).
        sim_sesID (str, optional): Replacement simulated session ID
            substituted for `subID` in the target file name. Defaults to None.
    """
    # Loop through each file-to-be-symlinked
    for file in file_subdir_dict.keys():
        # Define sub-directory of file
        subdir = file_subdir_dict[file]
        # Define full path of source file
        src_file = os.path.join(base_src_dir, subdir, file)
        if bld and new_bld:
            subdir = subdir.replace(bld, new_bld)
        # Create target dir
        tgt_dir = os.path.join(base_tgt_dir, subdir)
        os.makedirs(tgt_dir, exist_ok=True)
        # Create symbolic link of target file
        # Swap `subID` for `sim_sesID`
        if subID and sim_sesID:
            file = file.replace(subID, sim_sesID)
        # Swap `bld` for `new_bld`
        if bld and new_bld:
            # Mechanism to prevent accidental bug
            # whereby accidentally replace `bld` in
            # subID name. For example, if `bld` = 001,
            # and `new_bld` = 002, want to avoid
            # bug whereby change subID from "sub-00123"
            # to "sub-00223".
            subID_match = re.search(r"(sub-\d+)", file)
            if subID_match:
                correct_subID = subID_match.group(1)
            file = file.replace(bld, new_bld)
            if subID_match:
                file = re.sub(r"sub-\d+", correct_subID, file)
        # Define target file full path
        tgt_file = os.path.join(tgt_dir, file)
        # Create symlink
        if os.path.exists(tgt_file) or os.path.islink(tgt_file):
            os.remove(tgt_file)
        os.symlink(src_file, tgt_file)
        # Check symbolic link of target file
        # successfully created
        assert os.path.islink(
            tgt_file), f"ERROR: Symlink {tgt_file} not created."
        assert os.path.exists(
            tgt_file
        ), f"ERROR: Symlink {tgt_file} broken (target file doesn't exist)."


def generate_symlink_sim_dir(timestamp, main_config_file, pl_config_file,
                             list_file, run_type, reduce_mem_flag):
    """Build a "CBIG preproc structured" symlink dir for each subID/rm_run_th/combo/sim_sesID.

    Creates a directory that contains symbolic links (rather than real file
    copies) of all the files needed downstream (QC outlier files, bold/surf
    NIFTIs, registration files, masks, and normalization/warp files), pointing
    from the directory storing the real simulation files and the original
    preprocessed bold dir. This avoids creating unnecessary duplicates of the
    real simulation files, since each `combo` could reuse the same underlying
    real simulation files under different repeats/combinations. Writes a
    `success.txt` flag file once symlinking for a subID's rm_run_th/combo
    completes.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV file listing rows to group by
            `dataset`, `subID`, `rm_run_th`, `combo`, `sim_sesID`, each
            describing one `sim_bld`/`bld` combination to symlink.
        run_type (str): Run type; must be "ngt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode;
            controls which set of bold files are symlinked.
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    func2T1_sm = pl_config['func2T1']['func2T1_sm']
    if run_type == "ngt":
        pass
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'ngt'.")

    # Define types
    preproc_type = "after_mc"

    # Read in dataframe
    df = pd.read_csv(list_file, dtype={'bld': str, 'sim_bld': str})

    # Loop through each `rm_run_th`, `combo`, & `sim_sesID`
    for identifier, group_df in df.groupby(
            ['dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID']):

        # Extract identifiers
        dataset, subID, rm_run_th, combo, sim_sesID = identifier[
            0], identifier[1], identifier[2], identifier[3], identifier[4]
        unique_sim_bld_list = sorted(group_df['sim_bld'].unique())

        # Extract `REG_stem`
        resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                           pl_config_file,
                                                           timestamp,
                                                           preproc_type,
                                                           dataset,
                                                           subID,
                                                           run_type=run_type)
        BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                      resume_config_file)
        REG_stem = util.generate_pre_mc_stem_from_precursor_stem(BOLD_stem)
        REG_stem = util.append_suffix_to_stem(REG_stem, suffix="_mc_reg")

        # Cast to string
        rm_run_th = str(rm_run_th)

        # Create dir storing symlinks
        base_tgt_dir = util.get_sim_ses_preproc_dir_symlink_path(
            main_config_file,
            pl_config_file,
            dataset,
            subID,
            rm_run_th=rm_run_th,
            combo=combo,
            sim_sesID=sim_sesID)
        os.makedirs(base_tgt_dir, exist_ok=True)

        # Create `.bold` file
        logs_dir = os.path.join(base_tgt_dir, "logs")
        os.makedirs(logs_dir, exist_ok=True)
        bold_file = os.path.join(logs_dir, f"{sim_sesID}.bold")
        with open(bold_file, "w") as f:
            f.write(" ".join(unique_sim_bld_list))
        # Checks `.bold` file successfully created to
        # avoid additional unnecessary `check_XXX` function.
        assert os.path.exists(
            bold_file), f"ERROR: File {bold_file} not created."

        # Loop across vary threshold-level dataframe
        for idx, row in group_df.iterrows():
            # Extract identifiers
            sim_run_bin = row["sim_run_bin"]
            sim_bld = row["sim_bld"]
            vary_FD_metric_th = row[f"vary_{FD_metric}_th"]
            vary_DV_metric_th = row[f"vary_{DV_metric}_th"]

            # CREATE SYMLINK FOR QC OUTLIERS FILES
            # Define source dir
            base_src_dir = util.get_sim_ses_preproc_dir_real_path(
                main_config_file, pl_config_file, dataset, subID, sim_run_bin,
                sim_sesID)
            file_subdir_dict = {
                f"{sim_sesID}_bld{sim_bld}{FD_stem}{vary_FD_metric_th}{DV_stem}{vary_DV_metric_th}_motion_outliers.txt":
                "qc"
            }
            generate_symlink(file_subdir_dict, base_src_dir, base_tgt_dir)

        # Loop across sim_bld-level dataframe
        for idx, row in group_df.iterrows():
            # Extract identifiers
            sim_run_bin = row["sim_run_bin"]
            sim_bld = row["sim_bld"]
            bld = row["bld"]

            # CREATE SYMLINK FOR FILES FROM REAL SIMULATION DIR
            # Define source dir
            base_src_dir = util.get_sim_ses_preproc_dir_real_path(
                main_config_file, pl_config_file, dataset, subID, sim_run_bin,
                sim_sesID)

            # Extract `BOLD_stem`
            resume_config_file = util.get_config_file_dir_path(
                main_config_file,
                pl_config_file,
                timestamp,
                preproc_type,
                dataset,
                subID,
                run_type=run_type)
            BOLD_stem, _ = util.generate_stem_from_config(
                pl_config_file, resume_config_file)
            BOLD_stem_pre_mc = util.generate_pre_mc_stem_from_precursor_stem(
                BOLD_stem)

            # Define files and respective subdirs
            if reduce_mem_flag == "1":
                file_subdir_dict = {
                    f"{sim_sesID}_bld{sim_bld}{BOLD_stem}_sm{func2T1_sm}_loosemask.nii.gz":
                    f"bold/{sim_bld}",
                    f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{FD_stem}":
                    "bold/mc",
                    f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{DV_stem}":
                    "bold/mc",
                    f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{DV_stem}":
                    "bold/mc",
                    f"lh.{sim_sesID}_bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz":
                    "surf",
                    f"rh.{sim_sesID}_bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz":
                    "surf"
                }
            elif reduce_mem_flag == "0":
                file_subdir_dict = {
                    f"{sim_sesID}_bld{sim_bld}{BOLD_stem}_T1_sm{func2T1_sm}_loosemask.nii.gz":
                    f"bold/{sim_bld}",
                    f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{FD_stem}":
                    "bold/mc",
                    f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{DV_stem}":
                    "bold/mc",
                    f"{sim_sesID}_bld{sim_bld}{BOLD_stem_pre_mc}_motion_outliers{DV_stem}":
                    "bold/mc",
                    f"lh.{sim_sesID}_bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz":
                    "surf",
                    f"rh.{sim_sesID}_bld{sim_bld}{BOLD_stem}_fs6_sm6.nii.gz":
                    "surf"
                }
            else:
                raise ValueError(
                    f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
                )

            # Generate symbolic links
            generate_symlink(file_subdir_dict, base_src_dir, base_tgt_dir)

            # CREATE SYMLINK FOR FILES FROM ORIGINAL BOLD DIR
            # Define source dir
            base_src_dir = util.get_preproc_dir_path(main_config_file,
                                                     pl_config_file,
                                                     preproc_type,
                                                     dataset,
                                                     subID,
                                                     run_type=run_type)

            # Define files and respective subdirs
            if reduce_mem_flag == "1":
                file_subdir_dict = {
                    f"gm_mask.bin.nii.gz": f"bold/{bld}/rest_mask",
                    f"gm_mask_{bld}_T2.bin.nii.gz": f"bold/{bld}/rest_mask",
                    f"loosebrainmask.bin.nii.gz": f"bold/{bld}/rest_mask",
                    f"loosebrainmask_{bld}_T2.bin.nii.gz":
                    f"bold/{bld}/rest_mask",
                    f"loosebrainmask_{bld}_T2.nii.gz": f"bold/{bld}/rest_mask",
                    f"{subID}_bld{bld}{REG_stem}.dat": f"bold/{bld}"
                }
            else:
                file_subdir_dict = {
                    f"gm_mask.bin.nii.gz": f"bold/{bld}/rest_mask",
                    f"loosebrainmask.bin.nii.gz": f"bold/{bld}/rest_mask",
                }

            # Generate symbolic links
            generate_symlink(file_subdir_dict,
                             base_src_dir,
                             base_tgt_dir,
                             bld=bld,
                             new_bld=sim_bld,
                             subID=subID,
                             sim_sesID=sim_sesID)

        # Loop across sim_sesID-level dataframe
        sim_sesID_df = group_df.drop_duplicates(
            subset=['dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID'])
        assert len(sim_sesID_df) == 1, (
            f"ERROR: expected 1 row for `sim_sesID_df` but got "
            f"{len(sim_sesID_df)} rows instead. `sim_sesID_df` "
            f"derived from `group_df` which should only contain 1 "
            f"sim_sesID.")

        for idx, row in sim_sesID_df.iterrows():
            # Extract identifiers
            sim_run_bin = row["sim_run_bin"]
            sim_bld = row["sim_bld"]
            bld = row["bld"]

            # CREATE SYMLINK FOR FILES FROM ORIGINAL BOLD DIR
            # Define source dir
            base_src_dir = util.get_preproc_dir_path(main_config_file,
                                                     pl_config_file,
                                                     preproc_type,
                                                     dataset,
                                                     subID,
                                                     run_type=run_type)

            # Define files and respective subdirs
            file_subdir_dict = {
                "norm.nii.gz": "vol",
                "norm_MNI1mm0GenericAffine.mat": "vol",
                "norm_MNI1mm1InverseWarp.nii.gz": "vol",
                "norm_MNI1mm1Warp.nii.gz": "vol"
            }

            # Generate symbolic links
            generate_symlink(file_subdir_dict, base_src_dir, base_tgt_dir)

    # Generate `finish.txt` to represent that generation
    # of symlink dir for this subID's rm_run_th's combo
    # is successful.
    subID_dir = util.get_sim_ses_preproc_dir_symlink_path(
        main_config_file, pl_config_file, dataset, subID)
    success_file_path = os.path.join(subID_dir, "success.txt")
    with open(success_file_path, "w") as f:
        pass


def submit_generate_symlink_sim_dir_job(timestamp, jobprefix, main_config_file,
                                        pl_config_file, dataset, subID,
                                        list_file, run_type, reduce_mem_flag):
    """Submit a PBS job that runs `generate_symlink_sim_dir` for one subject.

    `list_file` should comprise only the simulation runs of a single selected
    subject. Cleans up (restarts) the subject's symlink base directory before
    submitting.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        jobprefix (str): Prefix used to build the PBS job name.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        list_file (str): Path to a CSV file listing this subject's simulation
            runs to symlink.
        run_type (str): Run type; must be "ngt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.
    """
    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if run_type == "ngt":
        rep_TMS_Kong2025_symlink_dir = pl_config['path'][
            'rep_TMS_Kong2025_ngt_symlink_dir']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'ngt'.")

    # Define types
    preproc_type = "after_mc"

    # Determine file paths and directories
    resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                       pl_config_file,
                                                       timestamp,
                                                       preproc_type,
                                                       dataset,
                                                       subID,
                                                       run_type=run_type)
    std_dir = os.path.join(rep_TMS_Kong2025_symlink_dir, timestamp, "job")
    os.makedirs(std_dir, exist_ok=True)

    # Initialize input vars
    BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                  resume_config_file)
    util.generate_pre_mc_stem_from_precursor_stem(BOLD_stem)

    # Clean up old output(s)
    # Delete `sim_bld` bold dir
    sim_ses_preproc_symlink_base_dir = util.get_sim_ses_preproc_dir_symlink_path(
        main_config_file, pl_config_file, dataset, subID)
    util.restart_dir(sim_ses_preproc_symlink_base_dir)

    # PBS job parameters
    jobname = f"{jobprefix}_{dataset}_{subID}"
    stderr_file = os.path.join(std_dir, f"{dataset}_{subID}.err")
    stdout_file = os.path.join(std_dir, f"{dataset}_{subID}.out")

    # Directly construct the command script for job submission
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; "
        f"sys.path.append(\\\"{TAN2026_TMSMOTION_DIR}/TMS/Kong2025\\\"); "
        f"import CBIG_TMSmotion_Kong2025 as TMS_Kong2025; "
        f"TMS_Kong2025.generate_symlink_sim_dir("
        f"\\\"{timestamp}\\\", \\\"{main_config_file}\\\", \\\"{pl_config_file}\\\", "
        f"\\\"{list_file}\\\", \\\"{run_type}\\\", \\\"{reduce_mem_flag}\\\");' "
        f"source deactivate")

    # Submit job
    very_high_memory = "1G"
    high_memory = "1G"
    intermediate_memory = "1G"
    normal_memory = "1G"
    if dataset == "Newbold" and (subID == "sub-cast2" or subID == "sub-cast3"):
        memory = very_high_memory
    elif dataset == "MSC" or (dataset == "Newbold" and subID == "sub-cast1"):
        memory = high_memory
    elif dataset == "Myconnectome":
        memory = intermediate_memory
    else:
        memory = normal_memory
    util.submit_job(walltime="3:00:00",
                    memory=memory,
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def check_generate_symlink_sim_dir_job(timestamp, try_count, main_config_file,
                                       pl_config_file, list_file):
    """Check whether a `generate_symlink_sim_dir` job completed successfully.

    Verifies that the `success.txt` flag file was written for the subID
    referenced in `list_file` (which is assumed to contain a single subID).

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number (accepted for interface consistency
            with other `check_XXX_job` functions; not currently used since
            success is checked via a flag file rather than per-attempt logs).
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV file listing this subject's simulation
            runs; must contain exactly one unique `dataset`/`subID`.

    Returns:
        tuple: (bool, str) whether the job succeeded, and `list_file` itself
        (returned as a stand-in for a failure list, since none is generated).
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rep_TMS_Kong2025_symlink_dir = pl_config['path'][
        'rep_TMS_Kong2025_ngt_symlink_dir']

    # Assume job successful
    job_success = True

    # Define input & output dir
    timestamp_dir = os.path.join(rep_TMS_Kong2025_symlink_dir, timestamp)
    os.makedirs(timestamp_dir, exist_ok=True)

    # Read in dataframe, which is assumed to
    # be on the subID-level already.
    df = pd.read_csv(list_file, dtype={'bld': str, 'sim_bld': str})

    # Loop through each subID
    for identifier, group_df in df.groupby(['dataset', 'subID']):
        # Extract identifiers
        dataset, subID = identifier[0], identifier[1]

        # Check if `success.txt` exists.
        subID_dir = util.get_sim_ses_preproc_dir_symlink_path(
            main_config_file, pl_config_file, dataset, subID)
        success_file_path = os.path.join(subID_dir, "success.txt")
        if not os.path.exists(success_file_path):
            job_success = False
        """
        # Check the log file for success
        success_text = "Exit Status: 0"
        job_success, msg = util.check_job_text(outSTATS_file, success_text)
        """
        # Break out of loop because `df` should
        # only have a single subID
        return job_success, list_file


#############################################################################
# ===========================================================================
# TMS_KONG2025 (BOTH NGT & GT)
# ===========================================================================
#############################################################################
def submit_TMS_Kong2025_job(timestamp,
                            jobprefix,
                            max_TMS_Kong2025_job,
                            main_config_file,
                            pl_config_file,
                            dataset,
                            subID,
                            run_type,
                            target_type,
                            rm_run_th=None,
                            combo=None,
                            sim_sesID=None,
                            vary_FD_th=None,
                            vary_DV_th=None,
                            reduce_mem_flag=None):
    """Submit the main TMS targeting MATLAB job (`CBIG_TMS_wrapper`) for one run.

    Resolves all the structural (reconall/betsurf), preprocessed, and TMS
    output directories needed for one run, cleans up any previous TMS output
    for it, and submits a MATLAB job that computes the TMS scalp target for
    that run. Behaves differently depending on `run_type`:

    - "ngt": operates on one simulated run, identified by `rm_run_th`,
      `combo`, `sim_sesID`, and the varied FD/DV thresholds; reads the
      preprocessed data from the symlinked simulation dir; passes an extra
      set of identifiers to `CBIG_TMS_wrapper` for output naming.
    - "gt": operates on one subject's single ground-truth run using the
      pipeline's fixed FD/DV thresholds; reads the preprocessed data from the
      split-bld preproc dir.

    Job submission is throttled to at most `max_TMS_Kong2025_job` concurrently
    queued/running jobs sharing `jobprefix` (enforced by
    `util.submit_TMS_Kong2025_job`).

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        jobprefix (str): Prefix used to build the PBS job name and to identify
            jobs sharing the same concurrency cap.
        max_TMS_Kong2025_job (int): Maximum number of concurrently
            queued/running jobs allowed for this `jobprefix`.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        run_type (str): Run type; "ngt" or "gt".
        target_type (str): TMS target/network type (e.g. "dep"), used to
            select memory allocation and namespace output directories.
        rm_run_th (str, optional): Run-removal threshold identifier (ngt only).
            Defaults to None.
        combo (str, optional): Combination identifier for the selected runs
            (ngt only). Defaults to None.
        sim_sesID (str, optional): Simulated session ID (ngt only). Defaults
            to None.
        vary_FD_th (float or str, optional): Varied FD threshold (ngt only).
            Defaults to None; for "gt" this is instead set internally from the
            pipeline's fixed `FD_th`.
        vary_DV_th (float or str, optional): Varied DV threshold (ngt only).
            Defaults to None; for "gt" this is instead set internally from the
            pipeline's fixed `DV_th`.
        reduce_mem_flag (str, optional): "1" for reduced-memory mode, "0" for
            full mode; affects the MATLAB command and memory allocation.
            Defaults to None.
    """

    # Read main config
    main_config = util.read_json_file(main_config_file)
    data_dir = main_config['path']['data_dir']
    reconall_map = main_config['preproc_struct']['reconall_map']
    first_N_T1_list = main_config['path']['first_N_T1_list']
    username = str(main_config['jobscheduler']['username'])

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    if run_type == "ngt":
        rep_TMS_Kong2025_TMS_dir = pl_config['path'][
            'rep_TMS_Kong2025_ngt_TMS_dir']
        FD_th = pl_config['preproc']['FD_no_interp_th']
        DV_th = pl_config['preproc']['DV_no_interp_th']
    elif run_type == "gt":
        rep_TMS_Kong2025_TMS_dir = pl_config['path'][
            'rep_TMS_Kong2025_gt_TMS_dir']
        FD_th = pl_config['preproc']['FD_th']
        DV_th = pl_config['preproc']['DV_th']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted values are 'ngt' or 'gt'."
        )

    # Define types
    preproc_type = "after_mc"
    rep_TMS_Kong2025_TMS_dir = os.path.join(rep_TMS_Kong2025_TMS_dir,
                                            target_type)
    TMS_type = f'Kong2025_{target_type}'

    # Structural-based dirs
    reconall_dataset, reconall_subID = util.extract_reconall_map(
        reconall_map, dataset, subID)
    reconall_sesID, reconall_runID = util.get_first_reconall_sesID_runID(
        reconall_dataset, reconall_subID, first_N_T1_list)
    reconall_base_dir = util.get_reconall_base_dir_path(
        reconall_dataset, reconall_subID, reconall_sesID, data_dir)
    reconall_dir = os.path.join(reconall_base_dir, reconall_runID)
    betsurf_dir = util.get_betsurf_dir_path(reconall_dataset, reconall_subID,
                                            reconall_sesID, reconall_runID,
                                            data_dir)

    # Get config file to extract `BOLD_stem` later
    resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                       pl_config_file,
                                                       timestamp,
                                                       preproc_type,
                                                       dataset,
                                                       subID,
                                                       run_type=run_type)

    # Adjust ground truth vs non-ground truth
    if run_type == "ngt":
        # Init job output dir
        std_dir = os.path.join(
            rep_TMS_Kong2025_TMS_dir, timestamp, "job", dataset, subID,
            rm_run_th, combo, sim_sesID,
            f"{FD_metric}{vary_FD_th}_{DV_metric}{vary_DV_th}")

        # Init preproc dir
        preproc_dir = util.get_sim_ses_preproc_dir_symlink_path(
            main_config_file,
            pl_config_file,
            dataset,
            subID,
            rm_run_th=rm_run_th,
            combo=combo,
            sim_sesID=sim_sesID)

        # Init TMS dir
        TMS_Kong2025_dir = util.get_TMS_dir_path(run_type,
                                                 main_config_file,
                                                 dataset,
                                                 pl_name,
                                                 subID,
                                                 rm_run_th,
                                                 combo,
                                                 sim_sesID,
                                                 FD_metric,
                                                 vary_FD_th,
                                                 DV_metric,
                                                 vary_DV_th,
                                                 TMS_type=TMS_type)

        # Init BOLD_stem
        BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                      resume_config_file)

        # Init PBS job parameters
        jobname = f"{jobprefix}_{dataset}_{subID}_{rm_run_th}_{combo}_{sim_sesID}_{vary_FD_th}_{vary_DV_th}"
        stderr_file = os.path.join(
            std_dir,
            f"{dataset}_{subID}_{rm_run_th}_{combo}_{sim_sesID}_{FD_metric}{vary_FD_th}_{DV_metric}{vary_DV_th}.err"
        )
        stdout_file = os.path.join(
            std_dir,
            f"{dataset}_{subID}_{rm_run_th}_{combo}_{sim_sesID}_{FD_metric}{vary_FD_th}_{DV_metric}{vary_DV_th}.out"
        )
    elif run_type == "gt":
        # Set previously used strictest thresholds as "varied" thresholds
        vary_FD_th = FD_th
        vary_DV_th = DV_th

        # Init job output dir
        std_dir = os.path.join(rep_TMS_Kong2025_TMS_dir, timestamp, "job",
                               dataset, subID)

        # Init preproc dir
        preproc_dir = util.get_split_bld_preproc_dir_path(
            main_config_file, pl_config_file, dataset, subID)

        # Init TMS dir
        TMS_Kong2025_dir = util.get_TMS_dir_path(run_type,
                                                 main_config_file,
                                                 dataset,
                                                 pl_name,
                                                 subID,
                                                 TMS_type=TMS_type)

        # Init BOLD_stem
        BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                      resume_config_file)

        # Init PBS job parameters
        jobname = f"{jobprefix}_{dataset}_{subID}"
        stderr_file = os.path.join(std_dir, f"{dataset}_{subID}.err")
        stdout_file = os.path.join(std_dir, f"{dataset}_{subID}.out")
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted values are 'ngt' or 'gt'."
        )

    # Create job output dir
    os.makedirs(std_dir, exist_ok=True)

    # Preproc-based dirs
    mask_dir = os.path.join(preproc_dir, "bold", "001", "rest_mask")
    assert os.path.exists(mask_dir), f"ERROR: {mask_dir} does not exist."

    # TMS-based output dirs
    out_roi_dir = os.path.join(TMS_Kong2025_dir, "ROI")
    out_fc_dir = os.path.join(TMS_Kong2025_dir, "FC")
    out_mshbm_dir = os.path.join(TMS_Kong2025_dir, "MSHBM")
    out_tms_dir = os.path.join(TMS_Kong2025_dir, "TMS")

    # Extract pre_FD & post_DV stems
    pre_FD_stem, post_DV_stem = util.generate_pre_FD_and_post_DV_stem(
        BOLD_stem, FD_stem, DV_stem)

    # Clean up old output(s)
    util.restart_dir(TMS_Kong2025_dir)

    # Directly construct the command script for job submission
    if run_type == "ngt":
        cmd_script = (
            f"matlab -nosplash -nodisplay -nodesktop -r \\\"clear; clc; close all; "
            f"addpath(\'{KONG2026_TMS_DIR}/targeting\'); "
            f"CBIG_TMS_wrapper(\'{reconall_dir}\', \'{preproc_dir}\', "
            f"\'{mask_dir}\', \'{betsurf_dir}\', \'{out_roi_dir}\', \'{out_fc_dir}\', "
            f"\'{out_mshbm_dir}\', \'{out_tms_dir}\', \'{FD_th}\',  \'{DV_th}\', "
            f"\'{pre_FD_stem}\', \'{FD_stem}\', \'{DV_stem}\',\'{post_DV_stem}\', "
            f"\'{vary_FD_th}\', \'{vary_DV_th}\', \'{pl_name}\', \'{run_type}\', "
            f"\'{reduce_mem_flag}\',\'{target_type}\',\'{rm_run_th}\', \'{combo}\', \'{sim_sesID}\'); "
            f"exit;\\\"")
    else:
        cmd_script = (
            f"matlab -nosplash -nodisplay -nodesktop -r \\\"clear; clc; close all; "
            f"addpath(\'{KONG2026_TMS_DIR}/targeting\'); "
            f"CBIG_TMS_wrapper(\'{reconall_dir}\', \'{preproc_dir}\', "
            f"\'{mask_dir}\', \'{betsurf_dir}\', \'{out_roi_dir}\', \'{out_fc_dir}\', "
            f"\'{out_mshbm_dir}\', \'{out_tms_dir}\', \'{FD_th}\',  \'{DV_th}\', "
            f"\'{pre_FD_stem}\', \'{FD_stem}\', \'{DV_stem}\',\'{post_DV_stem}\', "
            f"\'{vary_FD_th}\', \'{vary_DV_th}\', \'{pl_name}\', \'{run_type}\', "
            f"\'{reduce_mem_flag}\',\'{target_type}\'); exit;\\\"")

    # Define memory & walltime
    if reduce_mem_flag == "1":
        mem_dict = {
            "YaleTestRetest": {
                "dep": "18G"
            },
            "Newbold": {
                "dep": "10G"
            },
            "MSC": {
                "dep": "10G"
            },
            "NSD": {
                "dep": "16G"
            },
            "Myconnectome": {
                "dep": "16G"
            },
            "Kirby": {
                "dep": "10G"
            },
            "IPCAS6": {
                "dep": "9G"
            },
            "HBNSSI": {
                "dep": "9G"
            }
        }
        default_memory = "16G"

        if run_type == "ngt":
            walltime = "2:00:00"
            if dataset in mem_dict:
                memory = mem_dict[dataset][target_type]
            else:
                memory = default_memory
        elif run_type == "gt":
            high_memory = "30G"
            intermediate_memory = "25G"
            normal_memory = "20G"
            walltime = "6:00:00"
            if (dataset == "NSD") or (dataset
                                      == "Newbold") or (dataset
                                                        == "YaleTestRetest"):
                memory = high_memory
            elif (dataset == "HBNSSI" and subID == "sub-0031125") or (
                    dataset == "Myconnectome"):
                memory = intermediate_memory
            else:
                memory = normal_memory
        else:
            raise ValueError(
                f"Invalid run_type '{run_type}'. Accepted values are 'ngt' or 'gt'."
            )
    elif reduce_mem_flag == "0":
        high_memory = "100G"
        normal_memory = "40G"
        memory = high_memory if dataset in ["Newbold"] else normal_memory
        if run_type == "ngt":
            walltime = "2:00:00"
        elif run_type == "gt":
            walltime = "6:00:00"
        else:
            raise ValueError(
                f"Invalid run_type '{run_type}'. Accepted values are 'ngt' or 'gt'."
            )
    else:
        raise ValueError(
            f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
        )

    # Submit job
    util.submit_TMS_Kong2025_job(walltime=walltime,
                                 memory=memory,
                                 ncpus="1",
                                 jobname=jobname,
                                 stderr_file=stderr_file,
                                 stdout_file=stdout_file,
                                 cmd_script=cmd_script,
                                 max_jobs=max_TMS_Kong2025_job,
                                 username=username,
                                 jobprefix=jobprefix)


def submit_all_TMS_Kong2025_job(timestamp, jobprefix, max_TMS_Kong2025_job,
                                main_config_file, pl_config_file, run_type,
                                list_file, target_type, reduce_mem_flag):
    """Submit TMS_Kong2025 jobs for every simulated run listed in `list_file` (ngt only).

    Deduplicates `list_file` down to one row per unique simulated session and
    FD/DV threshold combination, then submits one job per row via
    `submit_TMS_Kong2025_job`. For "gt" run_type, exits the script (ground
    truth runs must instead call `submit_TMS_Kong2025_job` directly, since
    each subject only has a single ground-truth job).

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        jobprefix (str): Prefix used to build PBS job names.
        max_TMS_Kong2025_job (int): Maximum number of concurrently
            queued/running jobs allowed for this `jobprefix`.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        run_type (str): Run type; must be "ngt" (otherwise the script exits).
        list_file (str): Path to a CSV file listing simulation runs and their
            varied FD/DV thresholds.
        target_type (str): TMS target/network type used to select memory
            allocation and namespace output directories.
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']

    if run_type == "ngt":
        # Read list of simulation runs across all thresholds across all subIDs
        sim_run_df = pd.read_csv(list_file,
                                 dtype={
                                     'rm_run_th': str,
                                     'bld': str,
                                     'sim_bld': str,
                                     'highest_sim_bld': str
                                 })

        # Within each simulation session, only keep 1
        # simulation run (per FD & DVthreshold).
        sim_run_df = sim_run_df.drop_duplicates(subset=[
            'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
            f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
        ])

        # Only loop through each simulation session
        # (per FD & DV threshold) once.
        for _, row in sim_run_df.iterrows():
            dataset, subID, rm_run_th, combo, sim_sesID, vary_FD_th, vary_DV_th = row[
                [
                    'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
                    f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
                ]]

            # Submit job for 1 simulation run for 1 varied threshold
            submit_TMS_Kong2025_job(timestamp, jobprefix, max_TMS_Kong2025_job,
                                    main_config_file, pl_config_file, dataset,
                                    subID, run_type, target_type, rm_run_th,
                                    combo, sim_sesID, vary_FD_th, vary_DV_th,
                                    reduce_mem_flag)
    else:
        print(
            "Ground truth case will use `submit_TMS_Kong2025_job` instead. Exiting script now..."
        )
        sys.exit(1)


def check_TMS_Kong2025_job(timestamp,
                           try_count,
                           main_config_file,
                           pl_config_file,
                           data_dir,
                           dataset,
                           subID,
                           run_type,
                           target_type,
                           all_job_success=None,
                           rm_run_th=None,
                           combo=None,
                           sim_sesID=None,
                           highest_new_bld=None,
                           vary_FD_th=None,
                           vary_DV_th=None,
                           reduce_mem_flag=None):
    """Check whether one TMS_Kong2025 job's output files were created successfully.

    Verifies that the MSHBM/TMS target result files (and, for `reduce_mem_flag`
    "0", the per-bld SGC FC NIFTI files) expected from `submit_TMS_Kong2025_job`
    exist. Any missing file is logged to a per-attempt failure message file,
    and for "ngt" the corresponding identifiers are also logged to a
    per-attempt failure list file.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number, used to namespace the failure log files.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        data_dir (str): Data directory (read but unused beyond being passed
            through).
        dataset (str): Dataset name.
        subID (str): Subject ID.
        run_type (str): Run type; "ngt" or "gt".
        target_type (str): TMS target/network type used to namespace output
            directories.
        all_job_success (bool, optional): Running success flag accumulated
            across calls (ngt only); ignored and reset to True internally for
            "gt", since each gt subject has only one job. Defaults to None.
        rm_run_th (str, optional): Run-removal threshold identifier (ngt only).
            Defaults to None.
        combo (str, optional): Combination identifier for the selected runs
            (ngt only). Defaults to None.
        sim_sesID (str, optional): Simulated session ID (ngt only). Defaults
            to None.
        highest_new_bld (str, optional): Highest simulated bold ID for this
            session, used to determine how many per-bld SGC FC files to check
            when `reduce_mem_flag` is "0". Defaults to None.
        vary_FD_th (float or str, optional): Varied FD threshold (ngt only).
            Defaults to None.
        vary_DV_th (float or str, optional): Varied DV threshold (ngt only).
            Defaults to None.
        reduce_mem_flag (str, optional): "1" for reduced-memory mode (skip
            per-bld SGC FC file checks), "0" for full mode (check them).
            Defaults to None.

    Returns:
        bool or tuple: For run_type "gt", the `all_job_success` bool. For
        run_type "ngt", a tuple (bool, str) of the updated `all_job_success`
        flag and the path to the failure list file for this `try_count`.
    """
    # If ground truth, define all_job_success, because not using `check_all_TMS_Kong2025_job`
    # function (which defines all_job_success & then pass in to `check_TMS_Kong2025_job`).
    # Not using `check_all_TMS_Kong2025` because no list to loop through, as each subID =
    # only 1 ground truth TMS job)
    if run_type == "gt":
        all_job_success = True

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    MSHBM_w = pl_config['MSHBM']['MSHBM_w']
    MSHBM_MRF = pl_config['MSHBM']['MSHBM_MRF']
    if run_type == "ngt":
        rep_TMS_Kong2025_TMS_dir = pl_config['path'][
            'rep_TMS_Kong2025_ngt_TMS_dir']
    elif run_type == "gt":
        rep_TMS_Kong2025_TMS_dir = pl_config['path'][
            'rep_TMS_Kong2025_gt_TMS_dir']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted values are 'ngt' or 'gt'."
        )

    # Define types
    TMS_type = f'Kong2025_{target_type}'

    # Define files to check
    if run_type == "ngt":
        dir_name = sim_sesID
        TMS_Kong2025_dir = util.get_TMS_dir_path(run_type,
                                                 main_config_file,
                                                 dataset,
                                                 pl_name,
                                                 subID,
                                                 rm_run_th,
                                                 combo,
                                                 sim_sesID,
                                                 FD_metric,
                                                 vary_FD_th,
                                                 DV_metric,
                                                 vary_DV_th,
                                                 TMS_type=TMS_type)
    elif run_type == "gt":
        dir_name = subID
        TMS_Kong2025_dir = util.get_TMS_dir_path(run_type,
                                                 main_config_file,
                                                 dataset,
                                                 pl_name,
                                                 subID,
                                                 TMS_type=TMS_type)
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted values are 'ngt' or 'gt'."
        )

    filename2check = [
        f"TMS/target_results/group_weight/{dir_name}/MSHBM_scalp_target_results/MSHBM_target_results.mat",
        f"TMS/target_results/group_weight/{dir_name}/MSHBM_scalp_target_results/optimal_target_results.mat",
        f"TMS/target_results/group_weight/{dir_name}/MSHBM_scalp_target_results/scalp_dist_target.mat",
        f"MSHBM/ind_parcellation/test_set/Ind_parcellation_MSHBM_sub1_w{MSHBM_w}_MRF{MSHBM_MRF}.mat"
    ]

    if reduce_mem_flag == "1":
        pass
    elif reduce_mem_flag == "0":
        for i in range(1, int(highest_new_bld) + 1):
            new_bld_str = str(i).zfill(3)
            filename2check.append(
                f"TMS/Ind_SGC_FC/group_weight/{dir_name}/SGC_FC_group_weight_{new_bld_str}.nii.gz"
            )
            filename2check.append(
                f"TMS/Ind_SGC_FC/group_weight/{dir_name}/SGC_FC_group_weight_{new_bld_str}_DLPFC.nii.gz"
            )
    else:
        raise ValueError(
            f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
        )

    # Define output files
    timestamp_dir = os.path.join(rep_TMS_Kong2025_TMS_dir, timestamp,
                                 f"{dataset}_{subID}")
    os.makedirs(timestamp_dir, exist_ok=True)
    fail_list = os.path.join(timestamp_dir, f"fail_list_try_{try_count}.txt")
    fail_msg = os.path.join(timestamp_dir, f"fail_msg_try_{try_count}.txt")

    # Adjust output file header & data format according to non-ground truth vs ground truth
    if run_type == "ngt":
        fail_msg_header = (
            f"dataset,subID,rm_run_th,combo,sim_sesID,highest_sim_bld,vary_{FD_metric}_th,"
            f"vary_{DV_metric}_th    error_message\n")
        fail_msg_data = (
            f"{dataset},{subID},{rm_run_th},{combo},{sim_sesID},{highest_new_bld},"
            f"{vary_FD_th},{vary_DV_th}  {{filename}}\n")
        fail_list_header = (
            f"dataset,subID,rm_run_th,combo,sim_sesID,highest_sim_bld,vary_{FD_metric}_th,"
            f"vary_{DV_metric}_th\n")
        fail_list_data = (
            f"{dataset},{subID},{rm_run_th},{combo},{sim_sesID},{highest_new_bld},"
            f"{vary_FD_th},{vary_DV_th}\n")
    elif run_type == "gt":
        fail_msg_header = f"dataset,subID    error_message\n"
        fail_msg_data = f"{dataset},{subID}    {{filename}}\n"
        fail_list_header = f"dataset,subID\n"
        fail_list_data = f"{dataset},{subID}\n"
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted values are 'ngt' or 'gt'."
        )

    # Check job output files for success
    id_success = True
    for filename in filename2check:
        exist = util.check_job_output(os.path.join(TMS_Kong2025_dir, filename))
        if not exist:
            id_success = False
            # Handle failure message list
            if not os.path.exists(fail_msg):
                with open(fail_msg, 'w') as f:
                    f.write(fail_msg_header)
            with open(fail_msg, 'a') as f:
                f.write(fail_msg_data.format(filename=filename))

    if not id_success:
        all_job_success = False
        # Handle failed ID list
        if not os.path.exists(fail_list):
            with open(fail_list, 'w') as f:
                f.write(fail_list_header)
        with open(fail_list, 'a') as f:
            f.write(fail_list_data)

    if run_type == "ngt":
        return all_job_success, fail_list
    elif run_type == "gt":
        return all_job_success
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted values are 'ngt' or 'gt'."
        )


def check_all_TMS_Kong2025_job(timestamp, try_count, main_config_file,
                               pl_config_file, data_dir, list_file, run_type,
                               target_type, reduce_mem_flag):
    """Check every TMS_Kong2025 job (per simulated session/threshold) in `list_file` (ngt only).

    Deduplicates `list_file` down to one row per unique simulated session and
    FD/DV threshold combination, then checks each one via
    `check_TMS_Kong2025_job`. For "gt" run_type, exits the script (ground
    truth runs must instead call `check_TMS_Kong2025_job` directly).

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number, used to namespace the failure log files.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        data_dir (str): Data directory, passed through to
            `check_TMS_Kong2025_job`.
        list_file (str): Path to a CSV file listing simulation runs and their
            varied FD/DV thresholds.
        run_type (str): Run type; must be "ngt" (otherwise the script exits).
        target_type (str): TMS target/network type used to namespace output
            directories.
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.

    Returns:
        tuple: (bool, str) whether all checked jobs succeeded, and the path to
        the failure list file for this `try_count`.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']

    # Read list of simulation runs across all subIDs
    sim_run_df = pd.read_csv(list_file,
                             dtype={
                                 'rm_run_th': str,
                                 'bld': str,
                                 'sim_bld': str,
                                 'highest_sim_bld': str
                             })

    # Init check job vars
    all_job_success = True
    fail_list = None

    if run_type == "ngt":
        # Within each simulation session, only keep 1
        # simulation run (per FD & DV threshold).
        sim_run_df = sim_run_df.drop_duplicates(subset=[
            'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
            'highest_sim_bld', f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
        ])

        # Only loop through each simulation session
        # (per FD & DV threshold) once.
        for idx, row in sim_run_df.iterrows():
            dataset, subID, rm_run_th, combo, sim_sesID, highest_sim_bld, vary_FD_th, vary_DV_th = row[
                [
                    'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
                    'highest_sim_bld', f"vary_{FD_metric}_th",
                    f"vary_{DV_metric}_th"
                ]]

            # Check job
            all_job_success, fail_list = check_TMS_Kong2025_job(
                timestamp,
                try_count,
                main_config_file,
                pl_config_file,
                data_dir,
                dataset,
                subID,
                run_type,
                target_type,
                all_job_success=all_job_success,
                rm_run_th=rm_run_th,
                combo=combo,
                sim_sesID=sim_sesID,
                highest_new_bld=highest_sim_bld,
                vary_FD_th=vary_FD_th,
                vary_DV_th=vary_DV_th,
                reduce_mem_flag=reduce_mem_flag)
    else:
        print(
            "Ground truth case will use `check_TMS_Kong2025_job` instead. Exiting script now..."
        )
        sys.exit(1)

    return all_job_success, fail_list


#############################################################################
# ===========================================================================
# GT_LIST_2 GENERATION (ONLY GT)
# ===========================================================================
#############################################################################
def append_highest_split_bld_column(df):
    """Add a 'highest_split_bld' column containing each group's max 'split_bld'.

    Args:
        df (pd.DataFrame): DataFrame containing at least the columns
            'dataset', 'subID', and 'split_bld'.

    Returns:
        pd.DataFrame: `df` with a new 'highest_split_bld' column added, giving
        the highest 'split_bld' value within each ('dataset', 'subID') group.
    """

    # Create 'highest_split_bld' column
    df["highest_split_bld"] = df.groupby(["dataset", "subID"
                                          ])["split_bld"].transform("max")

    return df


def generate_gt_list_2_wrapper(timestamp, main_config_file, pl_config_file,
                               dataset, subID, run_type, reduce_mem_flag):
    """Process ground-truth data into `gt_list_2`, splitting long `bld`s into `split_bld`s.

    For a single selected subID, reads `gt_list_1` and processes each of the
    subject's ground-truth `bld`s: if the post-censored scan time of a `bld`
    exceeds `max_postcensored_min`, the `bld` is divided into multiple
    `split_bld`s with approximately evenly distributed non-censored frames;
    otherwise it becomes a single unsplit `split_bld`. Also generates
    symbolic links pointing from the original `bld` preproc folders (and the
    subject-level `vol` folder) to the new `split_bld` preproc folders.

    Note: although this wrapper only processes a single subID per call, it is
    named `generate_gt_list_2_wrapper` (rather than per-subID) because a full
    pipeline run invokes it once per selected subID to build up the complete
    `gt_list_2`.

    Side effects:
        - Saves the processed per-`split_bld` data for this subject as
          `gt_list_2.txt` in the output directory.
        - Writes a `.bold` file listing the subject's `split_bld`s and a
          `postcensored_scan_time.txt` file with total scan time.
        - Generates symbolic links for the required volume, rest-mask, and
          registration files in the output preprocessing directories.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        run_type (str): Run type; must be "gt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode;
            affects which rest-mask files are symlinked.
    """

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    rep_sel_subID_gt_dir = pl_config['path']['rep_sel_subID_gt_dir']
    rep_TMS_Kong2025_gt_gt_run_list_2_dir = pl_config['path'][
        'rep_TMS_Kong2025_gt_gt_run_list_2_dir']
    max_postcensored_min = pl_config['split']['max_postcensored_min']
    if run_type == "gt":
        FD_th = pl_config['preproc']['FD_th']
        DV_th = pl_config['preproc']['DV_th']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'gt'.")

    # Define types
    preproc_type = "after_mc"

    # Extract input lists
    gt_list_1 = os.path.join(rep_sel_subID_gt_dir, "gt_list_1.txt")

    # Convert into seconds
    max_postcensored_sec = max_postcensored_min * 60  # to seconds

    # Output dirs
    out_preproc_base_dir = util.get_split_bld_preproc_dir_path(
        main_config_file, pl_config_file, dataset, subID)
    out_logs_dir = os.path.join(out_preproc_base_dir, 'logs')
    if not os.path.exists(out_logs_dir):
        os.makedirs(out_logs_dir)

    # Read `bld`-level list
    df = pd.read_csv(gt_list_1, dtype={'bld': str, 'split_bld': str})
    df_subID = df[(df['dataset'] == dataset) & (df['subID'] == subID)]

    # Init vars for each subID, before each `bld`-level loop
    split_bld_value = 1
    df_split_blds = pd.DataFrame()
    split_blds = []
    vol_flag = 0

    # Create vol subfolder dict for each subID
    vol_dict = {
        "norm.nii.gz": "vol",
        "norm_MNI1mm1Warp.nii.gz": "vol",
        "norm_MNI1mm1InverseWarp.nii.gz": "vol",
        "norm_MNI1mm0GenericAffine.mat": "vol",
    }

    # Loop through each `bld`
    for runID_tuple, df_runIDs in df_subID.groupby(['dataset', 'subID',
                                                    'bld']):
        # Identifier variables
        bld = runID_tuple[2]

        # Check only have 1 original natural ground truth run
        assert len(
            df_runIDs
        ) == 1, f'ERROR: for {dataset}, {subID}, {bld}, > 1 unique `bld`!'

        # Define convenience variables
        TR = df_runIDs['TR'].iloc[0]

        # Define input dirs
        in_preproc_base_dir = util.get_preproc_dir_path(main_config_file,
                                                        pl_config_file,
                                                        preproc_type,
                                                        dataset,
                                                        subID,
                                                        run_type=run_type)

        # Define input files
        in_runID_outliers = os.path.join(
            in_preproc_base_dir, 'qc',
            f'{subID}_bld{bld}{FD_stem}{FD_th}{DV_stem}{DV_th}_motion_outliers.txt'
        )

        # Generate symbolic links for vol subdir
        if vol_flag == 0:
            generate_symlink(vol_dict, in_preproc_base_dir,
                             out_preproc_base_dir)
            vol_flag = 1
        elif vol_flag == 1:
            pass
        else:
            print(
                f'ERROR: for {dataset}, {subID}, unexpected value for vol_flag!'
            )
            sys.exit(1)

        # Read motion outliers file
        arr_runID_outliers = np.loadtxt(in_runID_outliers, dtype=int)
        arr_runID_outliers = np.reshape(arr_runID_outliers, (-1, 1))

        # Compute number of post-censored non-censored frames
        gt_nframes_postcensored_computed = np.sum(arr_runID_outliers == 1)

        # Compute post-censored ground truth scan time
        gt_run_time_postcensored_computed = gt_nframes_postcensored_computed * TR

        # If natural run time longer than threshold, split
        if gt_run_time_postcensored_computed > max_postcensored_sec:
            # Approximately evenly distribute non-censored frames across `split_bld`s
            gt_nframes_postcensored_th = max_postcensored_sec / TR
            num_split_blds = int(
                np.ceil(gt_nframes_postcensored_computed /
                        gt_nframes_postcensored_th))
            base_nframes = int(gt_nframes_postcensored_computed //
                               num_split_blds)
            remainder_nframes = int(gt_nframes_postcensored_computed %
                                    num_split_blds)
            split_blds_nframes = [base_nframes + 1] * remainder_nframes + [
                base_nframes
            ] * (num_split_blds - remainder_nframes)
            # Check approximately equal distribution of frames across `split_bld`s
            assert max(split_blds_nframes) - min(
                split_blds_nframes
            ) <= 1, f'for {dataset}, {subID}, {bld}, unequal distribution of frames across `split_bld`s!'
            # Check have exactly num_split_blds
            assert len(
                split_blds_nframes
            ) == num_split_blds, f'for {dataset}, {subID}, {bld}  have != num_split_blds!'
            # Check sum of number of frames of each split_bld tallies to exactly `bld` number of frames
            assert sum(
                split_blds_nframes
            ) == gt_nframes_postcensored_computed, f'for {dataset}, {subID}, {bld} missing frames in `split_bld`s!'
            # Record splitted
            split = 1
        # Else, don't split
        else:
            num_split_blds = 1
            split_blds_nframes = [gt_nframes_postcensored_computed]
            # Record not splitted
            split = 0

        # Record `split_bld`
        keep_indices = np.where(arr_runID_outliers == 1)[0]
        split_bld_position_start = 0  # inclusive
        split_bld_start_idx = keep_indices[
            split_bld_position_start]  # inclusive

        # Extract `REG_stem` for `.dat` file in `bld_level_dict`
        resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                           pl_config_file,
                                                           timestamp,
                                                           preproc_type,
                                                           dataset,
                                                           subID,
                                                           run_type=run_type)
        BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                      resume_config_file)
        REG_stem = util.generate_pre_mc_stem_from_precursor_stem(BOLD_stem)
        REG_stem = util.append_suffix_to_stem(REG_stem, suffix="_mc_reg")

        # Define bld-level files
        if reduce_mem_flag == "1":
            bld_level_dict = {
                f"gm_mask.bin.nii.gz": f"bold/{bld}/rest_mask",
                f"gm_mask_{bld}_T2.bin.nii.gz": f"bold/{bld}/rest_mask",
                f"loosebrainmask.bin.nii.gz": f"bold/{bld}/rest_mask",
                f"loosebrainmask_{bld}_T2.bin.nii.gz": f"bold/{bld}/rest_mask",
                f"loosebrainmask_{bld}_T2.nii.gz": f"bold/{bld}/rest_mask",
                f"{subID}_bld{bld}{REG_stem}.dat": f"bold/{bld}"
            }
        elif reduce_mem_flag == "0":
            bld_level_dict = {
                f"gm_mask.bin.nii.gz": f"bold/{bld}/rest_mask",
                f"loosebrainmask.bin.nii.gz": f"bold/{bld}/rest_mask",
                f"{subID}_bld{bld}{REG_stem}.dat": f"bold/{bld}"
            }
        else:
            raise ValueError(
                f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
            )

        for split_bld_nframes in split_blds_nframes:
            # Assign `split_bld`
            assert type(
                split_bld_value
            ) == int, f'for {dataset}, {subID}, {bld}, `split_bld_value`, {split_bld_value} NOT integer!'
            if split_bld_value < 10:
                split_bld = f'00{split_bld_value}'
            elif split_bld_value >= 10 and split_bld_value < 100:
                split_bld = f'0{split_bld_value}'
            else:
                split_bld = f'{split_bld_value}'
            split_blds.append(split_bld)

            # Generate symlink for rest_mask subdir
            generate_symlink(bld_level_dict,
                             in_preproc_base_dir,
                             out_preproc_base_dir,
                             bld=bld,
                             new_bld=split_bld)

            # Update split_bld
            assert type(
                split_bld_value
            ) == int, f'for {dataset}, {subID}, {bld}, `split_bld_value`, {split_bld_value} NOT integer!'
            split_bld_value += 1

            # Record split_bld_start_idx & split_bld_end_idx
            split_bld_position_end = split_bld_position_start + split_bld_nframes - 1  # inclusive
            split_bld_end_idx = keep_indices[
                split_bld_position_end]  # inclusive

            # Record gt_run_time_postcensored & gt_nframes_postcensored
            gt_nframes_postcensored = int(split_bld_position_end -
                                          split_bld_position_start + 1)
            gt_run_time_postcensored = round(gt_nframes_postcensored * TR, 1)

            # Add columns
            df_split_bld = df_runIDs.copy()
            new_cols = {
                'split_bld': split_bld,
                'split_bld_start_idx': split_bld_start_idx,
                'split_bld_end_idx': split_bld_end_idx,
                'split': split
            }
            for col, value in new_cols.items():
                assert col not in df_split_bld.columns, (
                    f'ERROR: for {dataset}, {subID}, {bld}, {split_bld} {col} column already exists!'
                )
                df_split_bld[col] = value

            cols2update = {
                'gt_run_time_postcensored': gt_run_time_postcensored,
                'gt_nframes_postcensored': gt_nframes_postcensored,
            }
            for col, value in cols2update.items():
                assert col in df_split_bld.columns, (
                    f'ERROR: for {dataset}, {subID}, {bld}, {split_bld} {col} column did not previously exist!'
                )
                df_split_bld[col] = value
            df_split_blds = pd.concat([df_split_blds, df_split_bld],
                                      ignore_index=True)

            # Update positions & idx
            # If have not reach end yet, keep sliding
            if split_bld_position_end < sum(split_blds_nframes) - 1:
                split_bld_position_start = split_bld_position_end + 1
                split_bld_start_idx = keep_indices[split_bld_position_start]
            # Elif reached end exactly
            elif split_bld_position_end == sum(split_blds_nframes) - 1:
                pass
                # print(f'{split_bld_position_start}_{split_bld_position_end}')
            # Unexpected case
            else:
                print(
                    f'ERROR: for {dataset}, {subID}, {bld}, unexpected case when determining '
                    f'`split_bld_start_idx` & `split_bld_end_idx`!')
                sys.exit(1)

    # Drop columns because not relevant to split_bld-level anymore
    col2drop = [
        f'gt_mean_{FD_metric}_postcensored',
        f'gt_mean_{DV_metric}_postcensored',
        f'gt_mean_{FD_metric}_znorm_postcensored',
        f'gt_mean_{DV_metric}_znorm_postcensored',
        'gt_mean_summary_postcensored'
    ]
    df_split_blds = util.drop_col(df_split_blds, col2drop)

    # Append `highest_split_bld` column
    df_split_blds = append_highest_split_bld_column(df_split_blds)

    # Reorder columns
    desired_columns = [
        'dataset', 'subID', 'split_bld', 'highest_split_bld', 'bld',
        'split_bld_start_idx', 'split_bld_end_idx', 'split', 'TR',
        'gt_run_time_postcensored', 'gt_nframes_postcensored'
    ]
    df_split_blds = util.reorder_df_columns(df_split_blds, desired_columns)

    # Check `split_bld` 1) ascending (increment by 1) & 2) unique values
    bld_values = df_split_blds['split_bld'].astype(int)

    # Check unique
    assert len(bld_values) == len(
        set(bld_values)), f'ERROR: for {dataset}, {subID}, {bld} \
        non-unique `split_bld` value present!'

    # Check that bld values MUST increment by 1
    # For example, should be bld 001 -> 002 -> 003 -> ...
    # For example, should not be bld 001 -> 005 -> ...
    if any(bld_values[i] != bld_values[i - 1] + 1
           for i in range(1, len(bld_values))):
        print(
            f'ERROR: for {dataset}, {subID}, {bld}, `split_bld` values NOT in ascending order or \
            do NOT increment by 1!')
        sys.exit(1)

    # Save df
    out_dir = os.path.join(rep_TMS_Kong2025_gt_gt_run_list_2_dir, "gt_list_2",
                           dataset, subID)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    df_split_blds.to_csv(os.path.join(out_dir, "gt_list_2.txt"), index=False)

    # Save .bold file
    with open(os.path.join(out_logs_dir, f'{subID}.bold'), 'w') as file:
        file.write(' '.join(split_blds))

    # Save total post-censored scan time
    total_hours = df_split_blds["gt_run_time_postcensored"].sum() / (60 * 60)
    output_file = os.path.join(out_dir, "postcensored_scan_time.txt")
    with open(output_file, "w") as f:
        f.write(f"{total_hours:.3f} hours\n")


def submit_generate_gt_list_2_job(timestamp, jobprefix, main_config_file,
                                  pl_config_file, dataset, subID, run_type,
                                  reduce_mem_flag):
    """Submit a PBS job that runs `generate_gt_list_2_wrapper` for one subject.

    Cleans up (restarts) the subject's split-bld preproc dir and gt_list_2
    output dir before submitting.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        jobprefix (str): Prefix used to build the PBS job name.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        run_type (str): Run type; must be "gt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.
    """

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rep_TMS_Kong2025_gt_gt_run_list_2_dir = pl_config['path'][
        'rep_TMS_Kong2025_gt_gt_run_list_2_dir']

    # Determine file paths and directories
    std_dir = os.path.join(rep_TMS_Kong2025_gt_gt_run_list_2_dir, timestamp,
                           "job")
    os.makedirs(std_dir, exist_ok=True)

    # Clean up old output(s)
    sub_run_preproc_dir = util.get_split_bld_preproc_dir_path(
        main_config_file, pl_config_file, dataset, subID)
    util.restart_dir(sub_run_preproc_dir)
    gt_list_2_dir = os.path.join(rep_TMS_Kong2025_gt_gt_run_list_2_dir,
                                 "gt_list_2", dataset, subID)
    util.restart_dir(gt_list_2_dir)

    # PBS job parameters
    jobname = f"{jobprefix}_{dataset}_{subID}"
    stderr_file = os.path.join(std_dir, f"{dataset}_{subID}.err")
    stdout_file = os.path.join(std_dir, f"{dataset}_{subID}.out")

    # Directly construct the command script for job submission
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; sys.path.append(\\\"{TAN2026_TMSMOTION_DIR}/TMS/Kong2025\\\"); "
        f"import CBIG_TMSmotion_Kong2025 as TMS_Kong2025; "
        f"TMS_Kong2025.generate_gt_list_2_wrapper(\\\"{timestamp}\\\", \
             \\\"{main_config_file}\\\", \\\"{pl_config_file}\\\", \\\"{dataset}\\\", \
             \\\"{subID}\\\", \\\"{run_type}\\\", \\\"{reduce_mem_flag}\\\")'; "
        f"source deactivate")

    # Submit job
    util.submit_job(walltime="1:00:00",
                    memory="1G",
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def check_generate_gt_list_2_job(timestamp, try_count, main_config_file,
                                 pl_config_file, dataset, subID, run_type,
                                 reduce_mem_flag):
    """Check whether a `generate_gt_list_2_wrapper` job's output files were created.

    Verifies that `gt_list_2.txt`, `postcensored_scan_time.txt`, the `.bold`
    log file, and the symlinked volume files exist. If `gt_list_2.txt` exists,
    also reads it to determine the expected per-`split_bld` rest-mask (and,
    for `reduce_mem_flag` "1", registration `.dat`) files and checks those too.
    Any missing file is logged to a per-attempt failure message file.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number, used to namespace the failure log file.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        run_type (str): Run type; must be "gt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode;
            affects which per-`split_bld` rest-mask files are expected.

    Returns:
        bool: Whether all expected output files exist.
    """
    # Assume job successful
    all_job_success = True

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if run_type == "gt":
        rep_TMS_Kong2025_gt_gt_run_list_2_dir = pl_config['path'][
            'rep_TMS_Kong2025_gt_gt_run_list_2_dir']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'gt'.")

    # Define types
    preproc_type = "after_mc"

    # Define filepath2check
    out_preproc_base_dir = util.get_split_bld_preproc_dir_path(
        main_config_file, pl_config_file, dataset, subID)
    gt_list_2_dir = os.path.join(rep_TMS_Kong2025_gt_gt_run_list_2_dir,
                                 "gt_list_2", dataset, subID)

    filepath2check = [
        f"{gt_list_2_dir}/gt_list_2.txt",
        f"{gt_list_2_dir}/postcensored_scan_time.txt",
        f"{out_preproc_base_dir}/logs/{subID}.bold",
        f"{out_preproc_base_dir}/vol/norm.nii.gz",
        f"{out_preproc_base_dir}/vol/norm_MNI1mm1Warp.nii.gz",
        f"{out_preproc_base_dir}/vol/norm_MNI1mm1InverseWarp.nii.gz",
        f"{out_preproc_base_dir}/vol/norm_MNI1mm0GenericAffine.mat"
    ]

    # Check existence of gt_list_2, & use it to append
    # `split_bld`-level files to filepath2check list.
    # (because previous lists do not contain `split_bld` column)
    gt_list_2 = os.path.join(gt_list_2_dir, "gt_list_2.txt")
    exist = util.check_job_output(gt_list_2)
    if exist:
        # Load gt_list_2 as a dataframe
        gt_list_2_df = pd.read_csv(gt_list_2,
                                   dtype={
                                       'bld': str,
                                       'split_bld': str
                                   })
        if "split_bld" not in gt_list_2_df.columns:
            raise ValueError(f"'split_bld' column not found in {gt_list_2}.")

        # Extract `REG_stem` for `.dat` file in `split_bld_level_dict`
        resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                           pl_config_file,
                                                           timestamp,
                                                           preproc_type,
                                                           dataset,
                                                           subID,
                                                           run_type=run_type)
        BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                      resume_config_file)
        REG_stem = util.generate_pre_mc_stem_from_precursor_stem(BOLD_stem)
        REG_stem = util.append_suffix_to_stem(REG_stem, suffix="_mc_reg")

        # Add file paths for each `split_bld` value
        if reduce_mem_flag == "1":
            for split_bld in gt_list_2_df["split_bld"]:
                split_bld_level_dict = {
                    f"gm_mask.bin.nii.gz": f"bold/{split_bld}/rest_mask",
                    f"gm_mask_{split_bld}_T2.bin.nii.gz":
                    f"bold/{split_bld}/rest_mask",
                    f"loosebrainmask.bin.nii.gz":
                    f"bold/{split_bld}/rest_mask",
                    f"loosebrainmask_{split_bld}_T2.bin.nii.gz":
                    f"bold/{split_bld}/rest_mask",
                    f"loosebrainmask_{split_bld}_T2.nii.gz":
                    f"bold/{split_bld}/rest_mask",
                    f"{subID}_bld{split_bld}{REG_stem}.dat":
                    f"bold/{split_bld}"
                }
                for filename, subdir in split_bld_level_dict.items():
                    filepath2check.append(
                        os.path.join(out_preproc_base_dir, subdir, filename))
        elif reduce_mem_flag == "0":
            for split_bld in gt_list_2_df["split_bld"]:
                split_bld_level_dict = {
                    f"gm_mask.bin.nii.gz": f"bold/{split_bld}/rest_mask",
                    f"loosebrainmask.bin.nii.gz":
                    f"bold/{split_bld}/rest_mask",
                }
                for filename, subdir in split_bld_level_dict.items():
                    filepath2check.append(
                        os.path.join(out_preproc_base_dir, subdir, filename))
        else:
            raise ValueError(
                f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
            )

    # Define output files
    timestamp_dir = os.path.join(rep_TMS_Kong2025_gt_gt_run_list_2_dir,
                                 timestamp, f"{dataset}_{subID}")
    os.makedirs(timestamp_dir, exist_ok=True)
    fail_msg = os.path.join(timestamp_dir, f"fail_msg_try_{try_count}.txt")

    # Check job output files for success
    for filepath in filepath2check:
        exist = util.check_job_output(filepath)
        if not exist:
            all_job_success = False
            # Handle failure message list
            if not os.path.exists(fail_msg):
                with open(fail_msg, 'w') as f:
                    f.write("dataset,subID    error_message\n")
            with open(fail_msg, 'a') as f:
                f.write(f"{dataset},{subID}    {filepath}\n")

    return all_job_success


#############################################################################
# ===========================================================================
# SPLIT GROUND TRUTH RUNS INTO SUB-RUNS (ONLY GT)
# ===========================================================================
#############################################################################


def generate_preproc_nii_path_for_gt_split(in_preproc_dir, dict_preproc,
                                           in_preproc_file_name,
                                           out_preproc_dir):
    """Resolve the input and output NIFTI file paths for one gt-split file entry.

    Given a nested dict describing how one input preprocessed file maps to
    one output (split) file, looks up the relevant sub-directory/file names
    and builds the corresponding full input and output paths, creating the
    output sub-directory if it does not already exist.

    Args:
        in_preproc_dir (str): Base directory of the input (unsplit) preprocessed
            data.
        dict_preproc (dict): Nested mapping of the form
            {in_file_name: {in_subdir_name: {out_file_name: [out_subdir_name]}}}
            describing how each input file maps to its output file/sub-dir.
        in_preproc_file_name (str): Key into `dict_preproc` identifying which
            input file's path to resolve.
        out_preproc_dir (str): Base directory of the output (split) preprocessed
            data.

    Returns:
        tuple: (str, str) the full input NIFTI file path and the full output
        NIFTI file path.
    """
    # Extract convenience variables
    in_preproc_subdir_name = next(iter(dict_preproc[in_preproc_file_name]))
    out_preproc_file_name = next(
        iter(dict_preproc[in_preproc_file_name][in_preproc_subdir_name]))
    out_preproc_subdir_name = next(
        iter(dict_preproc[in_preproc_file_name][in_preproc_subdir_name]
             [out_preproc_file_name]))

    # Define input dirs
    in_preproc_subdir = os.path.join(in_preproc_dir, in_preproc_subdir_name)
    # Define input files
    in_preproc_nii = os.path.join(in_preproc_subdir, in_preproc_file_name)
    # Define output dirs
    out_preproc_subdir = os.path.join(out_preproc_dir, out_preproc_subdir_name)
    if not os.path.exists(out_preproc_subdir):
        os.makedirs(out_preproc_subdir)
    # Define output files
    out_preproc_nii = os.path.join(out_preproc_subdir, out_preproc_file_name)

    return in_preproc_nii, out_preproc_nii


def gt_split_wrapper(main_config_file, pl_config_file, dataset, subID,
                     split_bld, bld, split_bld_start_idx, split_bld_end_idx,
                     split, timestamp, run_type, reduce_mem_flag):
    """Create one `split_bld`'s preprocessed files from its parent ground-truth `bld`.

    Top-level wrapper (intended to be run inside a submitted PBS job via
    `submit_gt_split_job`) that produces the surface, bold, and motion-outlier
    files for one `split_bld` derived from a ground-truth run `bld`:

    - If `split` == "1" (this `bld` was divided into multiple `split_bld`s),
      extracts the frame range [split_bld_start_idx, split_bld_end_idx]
      (inclusive) from each NIFTI via `fslroi`, and slices the corresponding
      rows out of the `bld`'s motion-outlier file.
    - If `split` == "0" (this `bld` was not divided), symlinks the whole
      original `bld` NIFTIs and motion-outlier file directly as the
      `split_bld`'s files (since output identical to input).

    Args:
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        split_bld (str): Output split-run ID being created.
        bld (str): Parent ground-truth run (bold) ID being split from.
        split_bld_start_idx (int or str): Inclusive starting frame index of
            this `split_bld` within `bld`; converted to int internally.
        split_bld_end_idx (int or str): Inclusive ending frame index of this
            `split_bld` within `bld`; converted to int internally.
        split (str): "1" if `bld` was divided into multiple `split_bld`s
            (extract via `fslroi`), "0" if not (symlink the whole run).
        timestamp (str): Timestamp identifying this pipeline run.
        run_type (str): Run type; must be "gt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.
    """

    # Convert `split_bld_start_idx` & `split_bld_end_idx` to integers
    split_bld_start_idx = int(split_bld_start_idx)
    split_bld_end_idx = int(split_bld_end_idx)

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    func2T1_sm = pl_config['func2T1']['func2T1_sm']
    if run_type == "gt":
        FD_th = pl_config['preproc']['FD_th']
        DV_th = pl_config['preproc']['DV_th']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'gt'.")

    # Define types
    preproc_type = "after_mc"

    # Extract BOLD_stem
    resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                       pl_config_file,
                                                       timestamp,
                                                       preproc_type,
                                                       dataset,
                                                       subID,
                                                       run_type=run_type)
    BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                  resume_config_file)

    # Define input dirs
    in_preproc_dir = util.get_preproc_dir_path(main_config_file,
                                               pl_config_file,
                                               preproc_type,
                                               dataset,
                                               subID,
                                               run_type=run_type)

    # Define input files
    in_runID_outliers = os.path.join(
        in_preproc_dir, "qc",
        f"{subID}_bld{bld}{FD_stem}{FD_th}{DV_stem}{DV_th}_motion_outliers.txt"
    )

    # Define output dirs
    out_preproc_dir = util.get_split_bld_preproc_dir_path(
        main_config_file, pl_config_file, dataset, subID)
    out_qc_dir = os.path.join(out_preproc_dir, 'qc')
    if not os.path.exists(out_qc_dir):
        os.makedirs(out_qc_dir)

    # Init NIFTI dicts; structure = {in_file:{in_dir:{out_file:{out_dir}}}}
    dict_preproc = {}
    for hemi in ['lh', 'rh']:
        dict_preproc[f'{hemi}.{subID}_bld{bld}{BOLD_stem}_fs6_sm6.nii.gz'] = {
            os.path.join('surf'): {
                f'{hemi}.{subID}_bld{split_bld}{BOLD_stem}_fs6_sm6.nii.gz':
                [os.path.join('surf')]
            }
        }
    if reduce_mem_flag == "1":
        dict_key = f'{subID}_bld{bld}{BOLD_stem}_sm{func2T1_sm}_loosemask.nii.gz'
        split_key = f'{subID}_bld{split_bld}{BOLD_stem}_sm{func2T1_sm}_loosemask.nii.gz'
        dict_preproc[dict_key] = {
            os.path.join('bold', bld): {
                split_key: [os.path.join('bold', f'{split_bld}')]
            }
        }
    elif reduce_mem_flag == "0":
        dict_key = f'{subID}_bld{bld}{BOLD_stem}_T1_sm{func2T1_sm}_loosemask.nii.gz'
        split_key = f'{subID}_bld{split_bld}{BOLD_stem}_T1_sm{func2T1_sm}_loosemask.nii.gz'
        dict_preproc[dict_key] = {
            os.path.join('bold', bld): {
                split_key: [os.path.join('bold', f'{split_bld}')]
            }
        }
    else:
        raise ValueError(
            f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
        )

    # If `split_bld` to be splitted,
    if split == "1":
        # Split NIFTIs
        for in_preproc_file_name in dict_preproc.keys():
            in_preproc_nii, out_preproc_nii = generate_preproc_nii_path_for_gt_split(
                in_preproc_dir, dict_preproc, in_preproc_file_name,
                out_preproc_dir)
            # Extract `split_bld` via fslroi
            util.run_command_line(
                f'fslroi {in_preproc_nii} {out_preproc_nii} {split_bld_start_idx} '
                f'{split_bld_end_idx - split_bld_start_idx + 1}')
        # Split motion outliers file
        # Input files
        arr_runID_outliers = np.loadtxt(in_runID_outliers, dtype=int)
        arr_runID_outliers = np.reshape(arr_runID_outliers, (-1, 1))
        arr_split_bld_outliers = arr_runID_outliers[
            split_bld_start_idx:split_bld_end_idx + 1]
        np.savetxt(os.path.join(
            out_qc_dir,
            f'{subID}_bld{split_bld}{FD_stem}{FD_th}{DV_stem}{DV_th}_motion_outliers.txt'
        ), arr_split_bld_outliers, fmt='%d')

    # Elif `split_bld` not to be splitted,
    elif split == "0":
        # Create NIFTI symbolic links
        for in_preproc_file_name in dict_preproc.keys():
            in_preproc_nii, out_preproc_nii = \
                generate_preproc_nii_path_for_gt_split(in_preproc_dir, dict_preproc, in_preproc_file_name,
                                                       out_preproc_dir)
            src = os.path.join(in_preproc_nii)
            target = os.path.join(out_preproc_nii)
            if os.path.exists(target):
                os.remove(target)
            os.symlink(src, target)
        # Create motion outliers symlink
        # Input files
        src = in_runID_outliers
        # Output files
        target = os.path.join(
            out_qc_dir,
            f'{subID}_bld{split_bld}{FD_stem}{FD_th}{DV_stem}{DV_th}_motion_outliers.txt'
        )
        if os.path.exists(target):
            os.remove(target)
        os.symlink(src, target)

    # Unexpected case
    else:
        print(
            f'ERROR: for {dataset}, {subID}, {split_bld}, unexpected case for \'split\' column!'
        )


def submit_gt_split_job(timestamp, jobprefix, main_config_file, pl_config_file,
                        dataset, subID, split_bld, bld, split_bld_start_idx,
                        split_bld_end_idx, split, run_type, reduce_mem_flag):
    """Submit a PBS job that creates one `split_bld` via `gt_split_wrapper`.

    Deletes any stale outputs from a previous attempt (surface files, the
    motion-outlier QC file, and the bold NIFTI) before submitting. Memory is
    chosen based on `dataset`/`subID`.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        jobprefix (str): Prefix used to build the PBS job name.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        split_bld (str): Output split-run ID being created.
        bld (str): Parent ground-truth run (bold) ID being split from.
        split_bld_start_idx (int or str): Inclusive starting frame index of
            this `split_bld` within `bld`.
        split_bld_end_idx (int or str): Inclusive ending frame index of this
            `split_bld` within `bld`.
        split (str): "1" if `bld` was divided into multiple `split_bld`s, "0"
            if not.
        run_type (str): Run type; must be "gt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.
    """

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    func2T1_sm = pl_config['func2T1']['func2T1_sm']
    if run_type == "gt":
        rep_TMS_Kong2025_gt_gt_split_dir = pl_config['path'][
            'rep_TMS_Kong2025_gt_gt_split_dir']
        FD_th = pl_config['preproc']['FD_th']
        DV_th = pl_config['preproc']['DV_th']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'gt'.")

    # Define types
    preproc_type = "after_mc"

    # Determine file paths and directories
    resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                       pl_config_file,
                                                       timestamp,
                                                       preproc_type,
                                                       dataset,
                                                       subID,
                                                       run_type=run_type)
    std_dir = os.path.join(rep_TMS_Kong2025_gt_gt_split_dir, timestamp, "job")
    os.makedirs(std_dir, exist_ok=True)

    # Init input vars
    BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                  resume_config_file)

    # Clean up old output(s)
    split_bld_preproc_dir = util.get_split_bld_preproc_dir_path(
        main_config_file, pl_config_file, dataset, subID)

    # surf files
    surf_dir = os.path.join(split_bld_preproc_dir, "surf")
    surf_lh_nifti = os.path.join(
        surf_dir, f"lh.{subID}_bld{split_bld}{BOLD_stem}_fs6_sm6.nii.gz")
    surf_rh_nifti = os.path.join(
        surf_dir, f"rh.{subID}_bld{split_bld}{BOLD_stem}_fs6_sm6.nii.gz")
    util.delete_file(surf_lh_nifti)
    util.delete_file(surf_rh_nifti)
    # qc file
    qc_dir = os.path.join(split_bld_preproc_dir, "qc")
    outlier_file = os.path.join(
        qc_dir,
        f"{subID}_bld{split_bld}{FD_stem}{FD_th}{DV_stem}{DV_th}_motion_outliers.txt"
    )
    util.delete_file(outlier_file)
    # bold file
    if reduce_mem_flag == "1":
        bold_nifti_file = os.path.join(
            split_bld_preproc_dir, "bold", split_bld,
            f"{subID}_bld{split_bld}{BOLD_stem}_sm{func2T1_sm}_loosemask.nii.gz"
        )
    elif reduce_mem_flag == "0":
        bold_nifti_file = os.path.join(
            split_bld_preproc_dir, "bold", split_bld,
            f"{subID}_bld{split_bld}{BOLD_stem}_T1_sm{func2T1_sm}_loosemask.nii.gz"
        )
    else:
        raise ValueError(
            f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
        )
    util.delete_file(bold_nifti_file)

    # PBS job parameters
    jobname = f"{jobprefix}_{dataset}_{subID}_{split_bld}"
    stderr_file = os.path.join(std_dir, f"{dataset}_{subID}_{split_bld}.err")
    stdout_file = os.path.join(std_dir, f"{dataset}_{subID}_{split_bld}.out")

    # Directly construct the command script for job submission
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; sys.path.append(\\\"{TAN2026_TMSMOTION_DIR}/TMS/Kong2025\\\"); "
        f"import CBIG_TMSmotion_Kong2025 as TMS_Kong2025; "
        f"TMS_Kong2025.gt_split_wrapper(\\\"{main_config_file}\\\", \
            \\\"{pl_config_file}\\\",  \\\"{dataset}\\\", \\\"{subID}\\\", \
             \\\"{split_bld}\\\", \\\"{bld}\\\", \\\"{split_bld_start_idx}\\\", \
              \\\"{split_bld_end_idx}\\\", \\\"{split}\\\", \\\"{timestamp}\\\", \
              \\\"{run_type}\\\", \\\"{reduce_mem_flag}\\\")'; "
        f"source deactivate")

    # Submit job
    very_high_memory = "8G"
    high_memory = "7G"
    normal_memory = "6G"
    if dataset == "Newbold" and (subID == "sub-cast2" or subID == "sub-cast3"):
        memory = very_high_memory
    elif dataset == "MSC" or (dataset == "Newbold" and subID == "sub-cast1"):
        memory = high_memory
    else:
        memory = normal_memory
    util.submit_job(walltime="24:00:00",
                    memory=memory,
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def submit_all_gt_split_job(timestamp, jobprefix, main_config_file,
                            pl_config_file, run_type, list_file,
                            reduce_mem_flag):
    """Submit gt_split jobs for every `split_bld` listed in `list_file`.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        jobprefix (str): Prefix used to build PBS job names.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        run_type (str): Run type; must be "gt".
        list_file (str): Path to a CSV file (`gt_list_2`) listing `split_bld`s
            to submit gt_split jobs for.
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.
    """
    # Read list of sub runs for a single selected subID
    split_bld_df = pd.read_csv(list_file, dtype={'bld': str, 'split_bld': str})

    # Loop through each sub run
    for idx, row in split_bld_df.iterrows():
        dataset, subID, split_bld, bld, split_bld_start_idx, split_bld_end_idx, split = row[
            [
                'dataset', 'subID', 'split_bld', 'bld', 'split_bld_start_idx',
                'split_bld_end_idx', 'split'
            ]]

        # Submit job for 1 sub run
        submit_gt_split_job(timestamp, jobprefix, main_config_file,
                            pl_config_file, dataset, subID, split_bld, bld,
                            split_bld_start_idx, split_bld_end_idx, split,
                            run_type, reduce_mem_flag)


def check_gt_split_job(timestamp, try_count, main_config_file, pl_config_file,
                       dataset, subID, split_bld, bld, split_bld_start_idx,
                       split_bld_end_idx, split, all_job_success, run_type,
                       reduce_mem_flag):
    """Check whether one gt_split job's output files were created successfully.

    Verifies that the surface, motion-outlier, and bold files expected from a
    `submit_gt_split_job` run exist. Any missing file is logged to a
    per-attempt failure message file, and the corresponding identifiers are
    logged to a per-attempt failure list file.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number, used to namespace the failure log files.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset name.
        subID (str): Subject ID.
        split_bld (str): Output split-run ID.
        bld (str): Parent ground-truth run (bold) ID (used only in failure log
            messages).
        split_bld_start_idx (int or str): Inclusive starting frame index
            (used only in failure log messages).
        split_bld_end_idx (int or str): Inclusive ending frame index (used
            only in failure log messages).
        split (str): "1" or "0" (used only in failure log messages).
        all_job_success (bool): Running success flag accumulated across calls;
            set to False if this `split_bld`'s files are missing.
        run_type (str): Run type; must be "gt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.

    Returns:
        tuple: (bool, str) the updated `all_job_success` flag, and the path to
        the failure list file for this `try_count`.
    """

    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_stem = pl_config['stem']['FD_stem']
    DV_stem = pl_config['stem']['DV_stem']
    func2T1_sm = pl_config['func2T1']['func2T1_sm']
    if run_type == "gt":
        rep_TMS_Kong2025_gt_gt_split_dir = pl_config['path'][
            'rep_TMS_Kong2025_gt_gt_split_dir']
        FD_th = pl_config['preproc']['FD_th']
        DV_th = pl_config['preproc']['DV_th']
    else:
        raise ValueError(
            f"Invalid run_type '{run_type}'. Accepted value is only 'gt'.")

    # Define types
    preproc_type = "after_mc"

    # Define files to check
    resume_config_file = util.get_config_file_dir_path(main_config_file,
                                                       pl_config_file,
                                                       timestamp,
                                                       preproc_type,
                                                       dataset,
                                                       subID,
                                                       run_type=run_type)
    BOLD_stem, _ = util.generate_stem_from_config(pl_config_file,
                                                  resume_config_file)

    split_bld_preproc_dir = util.get_split_bld_preproc_dir_path(
        main_config_file, pl_config_file, dataset, subID)
    filepath2check = [
        f"{split_bld_preproc_dir}/surf/lh.{subID}_bld{split_bld}{BOLD_stem}_fs6_sm6.nii.gz",
        f"{split_bld_preproc_dir}/surf/rh.{subID}_bld{split_bld}{BOLD_stem}_fs6_sm6.nii.gz",
        f"{split_bld_preproc_dir}/qc/{subID}_bld{split_bld}{FD_stem}{FD_th}{DV_stem}{DV_th}_motion_outliers.txt"
    ]

    if reduce_mem_flag == "1":
        filepath2check.append(
            f"{split_bld_preproc_dir}/bold/{split_bld}/{subID}_bld{split_bld}{BOLD_stem}"
            f"_sm{func2T1_sm}_loosemask.nii.gz")
    elif reduce_mem_flag == "0":
        filepath2check.append(
            f"{split_bld_preproc_dir}/bold/{split_bld}/{subID}_bld{split_bld}{BOLD_stem}"
            f"_T1_sm{func2T1_sm}_loosemask.nii.gz")
    else:
        raise ValueError(
            f"Invalid value for reduce_mem_flag: {reduce_mem_flag}. Accepted values are '1' or '0'."
        )

    # Define output files
    timestamp_dir = os.path.join(rep_TMS_Kong2025_gt_gt_split_dir, timestamp,
                                 f"{dataset}_{subID}")
    os.makedirs(timestamp_dir, exist_ok=True)
    fail_list = os.path.join(timestamp_dir, f"fail_list_try_{try_count}.txt")
    fail_msg = os.path.join(timestamp_dir, f"fail_msg_try_{try_count}.txt")

    # Check job output files for success
    id_success = True
    for filepath in filepath2check:
        exist = util.check_job_output(os.path.join(filepath))
        if not exist:
            id_success = False
            # Handle failure message list
            if not os.path.exists(fail_msg):
                with open(fail_msg, 'w') as f:
                    f.write(
                        "dataset,subID,split_bld,bld,split_bld_start_idx,split_bld_end_idx,split    error_message\n"
                    )
            with open(fail_msg, 'a') as f:
                f.write(
                    f"{dataset},{subID},{split_bld},{bld},{split_bld_start_idx},"
                    f"{split_bld_end_idx},{split}    {filepath}\n")
    if not id_success:
        all_job_success = False
        # Handle failed ID list
        if not os.path.exists(fail_list):
            with open(fail_list, 'w') as f:
                f.write(
                    "dataset,subID,split_bld,bld,split_bld_start_idx,split_bld_end_idx,split\n"
                )
        with open(fail_list, 'a') as f:
            f.write(
                f"{dataset},{subID},{split_bld},{bld},{split_bld_start_idx},{split_bld_end_idx},{split}\n"
            )

    return all_job_success, fail_list


def check_all_gt_split_job(timestamp, try_count, main_config_file,
                           pl_config_file, list_file, run_type,
                           reduce_mem_flag):
    """Check every gt_split job (per `split_bld`) listed in `list_file`.

    Args:
        timestamp (str): Timestamp identifying this pipeline run.
        try_count (int): Attempt number, used to namespace the failure log files.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV file (`gt_list_2`) listing `split_bld`s
            to check.
        run_type (str): Run type; must be "gt".
        reduce_mem_flag (str): "1" for reduced-memory mode, "0" for full mode.

    Returns:
        tuple: (bool, str) whether all checked jobs succeeded, and the path to
        the failure list file for this `try_count` (None if `list_file` is empty).
    """

    # Read list of simulation runs across all subIDs
    split_bld_df = pd.read_csv(list_file, dtype={'bld': str, 'split_bld': str})

    # Init check job vars
    all_job_success = True
    fail_list = None

    # Loop through each simulation run
    for idx, row in split_bld_df.iterrows():
        dataset, subID, split_bld, bld, split_bld_start_idx, split_bld_end_idx, split = row[
            [
                'dataset', 'subID', 'split_bld', 'bld', 'split_bld_start_idx',
                'split_bld_end_idx', 'split'
            ]]

        # Check job
        all_job_success, fail_list = check_gt_split_job(
            timestamp, try_count, main_config_file, pl_config_file, dataset,
            subID, split_bld, bld, split_bld_start_idx, split_bld_end_idx,
            split, all_job_success, run_type, reduce_mem_flag)

    return all_job_success, fail_list
