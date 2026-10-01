'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

Examples:

python $TAN2026_TMSMOTION_DIR/compare/CBIG_TMSmotion_compare_Kong2025.py \
    $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json

python $TAN2026_TMSMOTION_DIR/compare/CBIG_TMSmotion_compare_Kong2025.py \
    $TAN2026_TMSMOTION_DIR/replication/pl2/config/pipeline.json

'''
#############################################################################
# ===========================================================================
# IMPORTS
# ===========================================================================
#############################################################################
import os
import pandas as pd
import sys
import numpy as np
import scipy.io

TAN2026_TMSMOTION_DIR = os.environ.get('TAN2026_TMSMOTION_DIR')
sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'util'))
import CBIG_TMSmotion_util as util

CBIG_CODE_DIR = util.get_cbig_code_dir()


#############################################################################
# ===========================================================================
# TMS TARGET EUCLIDEAN DISTANCE
# ===========================================================================
#############################################################################
def compare_Kong2025_TMS_mat2npy(target_mat_dir,
                                 target_npy_dir,
                                 TMS_type='Kong2025_dep'):
    """Convert a ground-truth or simulated TMS target .mat file to .npy.

    Reads `optimal_target_results.mat` from `target_mat_dir`, extracts the
    optimal TMS target centroid RAS coordinates, and saves them as
    `optimal_target_results.npy` in `target_npy_dir`. If the input .mat file
    is missing or cannot be parsed, a warning is printed and no .npy file is
    written.

    Args:
        target_mat_dir (str): Directory containing the input
            `optimal_target_results.mat` file.
        target_npy_dir (str): Directory to save the output
            `optimal_target_results.npy` file. Created if it does not exist.
        TMS_type (str, optional): TMS type identifier. Defaults to
            'Kong2025_dep', currently the only accepted value.
    """
    # Define paths
    target_mat_file = os.path.join(target_mat_dir,
                                   "optimal_target_results.mat")
    target_npy_file = os.path.join(target_npy_dir,
                                   "optimal_target_results.npy")
    os.makedirs(target_npy_dir, exist_ok=True)

    if os.path.exists(target_mat_file):
        # Read TMS target .mat
        target_mat = scipy.io.loadmat(target_mat_file)

        # Extract TMS target
        try:
            if TMS_type == 'Kong2025_dep':
                target_arr = target_mat['optimal_target_centroid_ras'][0]
            else:
                raise ValueError(
                    f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
                )

            # Save extracted TMS target
            np.save(target_npy_file, target_arr)
            # print(f"Done mat2npy for {target_npy_file}")
        except Exception as e:
            print(f"WARNING: {target_mat_file} problematic! Error: {e}")
    else:
        print(f"WARNING: {target_mat_file} does NOT exist!")


def compare_Kong2025_TMS_compute_euclidean_distance(RAS_1_npy, RAS_2_npy):
    """Compute the Euclidean distance between two TMS target coordinates.

    Args:
        RAS_1_npy (np.ndarray): First RAS coordinate array.
        RAS_2_npy (np.ndarray): Second RAS coordinate array.

    Returns:
        float: Euclidean (L2 norm) distance between `RAS_1_npy` and
        `RAS_2_npy`.
    """
    return np.linalg.norm(RAS_1_npy - RAS_2_npy)


def compare_Kong2025_TMS_wrapper(main_config_file,
                                 pl_config_file,
                                 list_file,
                                 mode,
                                 TMS_type='Kong2025_dep'):
    """Compute TMS target Euclidean distances for all subjects and sim runs.

    For each subject in `list_file`, converts the ground-truth TMS target
    .mat file to .npy and loads the ground-truth coordinates. Then, for
    every simulation run listed in that subject's `sim_run_list_3_{mode}.txt`
    (deduplicated to one run per simulation session and FD/DV threshold
    combination), converts the simulated TMS target .mat file to .npy,
    computes the Euclidean distance to the ground-truth coordinates, and
    saves it as `eucli_dist.npy` alongside the simulated .npy output.

    Args:
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV listing the selected dataset/subID
            combinations to process.
        mode (str): Simulation selection mode (e.g. 'grid' or 'slide') used
            to locate the `sim_run_list_3_{mode}.txt` file.
        TMS_type (str, optional): TMS type identifier. Defaults to
            'Kong2025_dep', currently the only accepted value.
    """
    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']

    ##########################################################################
    # Ground truth mat2npy
    ##########################################################################
    # Init output_type
    output_type = "TMS"

    # Loop through selected subID list
    selected_subID_df = pd.read_csv(list_file,
                                    dtype={
                                        'rm_run_th': str,
                                        'bld': str,
                                        'sim_bld': str,
                                        'highest_sim_bld': str,
                                        'split_bld': str
                                    })
    for idx, row in selected_subID_df.iterrows():
        # Init run_type
        run_type = "gt"

        # Define identifiers
        dataset, subID = row[['dataset', 'subID']]

        # Define input dirs
        TMS_Kong2025_dir = util.get_TMS_dir_path(run_type,
                                                 main_config_file,
                                                 dataset,
                                                 pl_name,
                                                 subID,
                                                 TMS_type=TMS_type)
        if TMS_type == 'Kong2025_dep':
            target_mat_dir = os.path.join(TMS_Kong2025_dir, "TMS",
                                          "target_results", "group_weight",
                                          subID, "MSHBM_scalp_target_results")
        else:
            raise ValueError(
                f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
            )

        # Define output dirs
        target_npy_dir = util.get_compare_Kong2025_output_dir_path(
            main_config_file,
            pl_config_file,
            output_type,
            dataset,
            subID,
            TMS_type=TMS_type)

        # Convert .mat to .npy files
        compare_Kong2025_TMS_mat2npy(target_mat_dir,
                                     target_npy_dir,
                                     TMS_type=TMS_type)

        # Load ground truth TMS coordinates
        RAS_gt_arr = np.load(
            os.path.join(target_npy_dir, "optimal_target_results.npy"))
        ##########################################################################
        # Simulation mat2npy
        ##########################################################################
        # Init run_type
        run_type = "ngt"

        # Loop through sim_run-level list
        sim_run_list_file = os.path.join(rep_sel_subID_sim_dir, "result",
                                         "summary", f"sim_run_list_3_{mode}",
                                         dataset, subID,
                                         f"sim_run_list_3_{mode}.txt")
        sim_run_df = pd.read_csv(sim_run_list_file,
                                 dtype={
                                     'rm_run_th': str,
                                     'bld': str,
                                     'sim_bld': str,
                                     'highest_sim_bld': str,
                                     'split_bld': str
                                 })

        # Within each simulation session, only keep 1
        # simulation run (per FD & DV threshold).
        sim_run_df = sim_run_df.drop_duplicates(subset=[
            'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
            f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
        ])

        # Loop through each simulation session
        for _, sim_run_row in sim_run_df.iterrows():
            # Define identifiers
            dataset, subID, rm_run_th, combo, sim_sesID, vary_FD_th, vary_DV_th = sim_run_row[
                [
                    'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
                    f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
                ]]

            # Define input dirs
            TMS_Kong2025_dir = util.get_TMS_dir_path(run_type,
                                                     main_config_file,
                                                     dataset,
                                                     pl_name,
                                                     subID,
                                                     rm_run_th=rm_run_th,
                                                     combo=combo,
                                                     sim_sesID=sim_sesID,
                                                     FD_metric=FD_metric,
                                                     vary_FD_th=vary_FD_th,
                                                     DV_metric=DV_metric,
                                                     vary_DV_th=vary_DV_th,
                                                     TMS_type=TMS_type)
            if TMS_type == 'Kong2025_dep':
                target_mat_dir = os.path.join(TMS_Kong2025_dir, "TMS",
                                              "target_results", "group_weight",
                                              sim_sesID,
                                              "MSHBM_scalp_target_results")
            else:
                raise ValueError(
                    f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
                )

            # Define output dirs
            target_npy_dir = util.get_compare_Kong2025_output_dir_path(
                main_config_file=main_config_file,
                pl_config_file=pl_config_file,
                output_type=output_type,
                dataset=dataset,
                subID=subID,
                rm_run_th=rm_run_th,
                combo=combo,
                sim_sesID=sim_sesID,
                FD_metric=FD_metric,
                vary_FD_th=vary_FD_th,
                DV_metric=DV_metric,
                vary_DV_th=vary_DV_th,
                TMS_type=TMS_type)

            # Convert .mat to .npy files
            compare_Kong2025_TMS_mat2npy(target_mat_dir,
                                         target_npy_dir,
                                         TMS_type=TMS_type)

            # Load Simulation TMS coordinates
            RAS_ngt_npy = os.path.join(target_npy_dir,
                                       "optimal_target_results.npy")

            # Compute and save Euclidean distance
            if os.path.exists(RAS_ngt_npy):
                RAS_ngt_arr = np.load(RAS_ngt_npy)
                eucli_dist = compare_Kong2025_TMS_compute_euclidean_distance(
                    RAS_gt_arr, RAS_ngt_arr)
                np.save(os.path.join(target_npy_dir, "eucli_dist.npy"),
                        eucli_dist)


def submit_compare_Kong2025_TMS_wrapper_job(timestamp, jobprefix,
                                            main_config_file, pl_config_file,
                                            list_file, mode, TMS_type):
    """Submit a PBS job that runs `compare_Kong2025_TMS_wrapper`.

    Args:
        timestamp (str): Timestamp string used to namespace the job's
            stdout/stderr directory.
        jobprefix (str): Prefix used to build the PBS job name.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV listing the selected dataset/subID
            combinations to process.
        mode (str): Simulation selection mode (e.g. 'grid' or 'slide').
        TMS_type (str): TMS type identifier, currently only 'Kong2025_dep'
            is a valid value.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )

    # Determine job std dir & flag file
    std_dir = os.path.join(rep_compare_Kong2025_dir, "TMS_wrapper", timestamp,
                           "job")
    os.makedirs(std_dir, exist_ok=True)
    flag_file = os.path.join(std_dir, "done.flag")

    # Construct command
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; sys.path.append(\\\"{TAN2026_TMSMOTION_DIR}/compare\\\"); "
        f"import CBIG_TMSmotion_compare_Kong2025 as compare_Kong2025; "
        f"compare_Kong2025.compare_Kong2025_TMS_wrapper(\\\"{main_config_file}\\\", \
            \\\"{pl_config_file}\\\", \\\"{list_file}\\\", \\\"{mode}\\\", \
            \\\"{TMS_type}\\\")'; "
        f"touch {flag_file}; "
        f"conda deactivate")

    # Init PBS job parameters
    jobname = f"{jobprefix}"
    stderr_file = os.path.join(std_dir, f"{jobprefix}.err")
    stdout_file = os.path.join(std_dir, f"{jobprefix}.out")

    # Submit job to headnode
    util.submit_job(walltime="24:00:00",
                    memory="4G",
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def check_compare_Kong2025_TMS_wrapper_job(timestamp, try_count,
                                           main_config_file, pl_config_file,
                                           TMS_type):
    """Check whether the `compare_Kong2025_TMS_wrapper` job completed.

    Args:
        timestamp (str): Timestamp string identifying the job run, used to
            locate its stdout/stderr directory.
        try_count (int): Retry attempt count. Not used directly in this
            function; kept for interface consistency with other
            `check_*_job` functions.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        TMS_type (str): TMS type identifier, currently only 'Kong2025_dep'
            is a valid value.

    Returns:
        bool: True if the job's flag file indicates successful completion,
        False otherwise.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )

    # Check flag file
    flag_file = os.path.join(rep_compare_Kong2025_dir, "TMS_wrapper",
                             timestamp, "job", "done.flag")
    return util.check_job_output(flag_file)


#############################################################################
# ===========================================================================
# MSHBM DICE COEFFICIENT
# ===========================================================================
#############################################################################
def compare_Kong2025_MSHBM_mat2npy(pl_config_file, label_mat_dir,
                                   label_npy_dir):
    """Convert an MSHBM individual parcellation .mat file to .npy files.

    Reads the MSHBM weighting (`MSHBM_w`) and MRF (`MSHBM_MRF`) parameters
    from the pipeline config to locate the input
    `Ind_parcellation_MSHBM_sub1_w{MSHBM_w}_MRF{MSHBM_MRF}.mat` file in
    `label_mat_dir`, extracts the left- and right-hemisphere label arrays,
    validates they are shaped (N, 1), flattens them to (N,), and saves them
    as separate .npy files in `label_npy_dir`. If the input .mat file is
    missing, a warning is printed and no .npy files are written.

    Args:
        pl_config_file (str): Path to the pipeline-specific config JSON
            file, used to look up `MSHBM_w` and `MSHBM_MRF`.
        label_mat_dir (str): Directory containing the input MSHBM
            parcellation .mat file.
        label_npy_dir (str): Directory to save the output lh/rh label .npy
            files. Created if it does not exist.
    """
    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    MSHBM_w = pl_config['MSHBM']['MSHBM_w']
    MSHBM_MRF = pl_config['MSHBM']['MSHBM_MRF']

    # Define paths
    os.makedirs(label_npy_dir, exist_ok=True)
    label_mat_file = os.path.join(
        label_mat_dir,
        f"Ind_parcellation_MSHBM_sub1_w{MSHBM_w}_MRF{MSHBM_MRF}.mat")
    lh_label_npy = os.path.join(label_npy_dir,
                                f"lh_labels_w{MSHBM_w}_MRF{MSHBM_MRF}.npy")
    rh_label_npy = os.path.join(label_npy_dir,
                                f"rh_labels_w{MSHBM_w}_MRF{MSHBM_MRF}.npy")

    if os.path.exists(label_mat_file):
        # Read MSHBM label .mat
        label_mat = scipy.io.loadmat(label_mat_file)

        # Extract left & right hemisphere labels
        lh_label_arr = label_mat['lh_labels']
        rh_label_arr = label_mat['rh_labels']

        # Assert lh_label_arr is of shape (N, 1)
        assert lh_label_arr.ndim == 2 and lh_label_arr.shape[1] == 1, \
            f"lh_label_arr must have shape (N, 1), but got {lh_label_arr.shape}"

        # Assert rh_label_arr is of shape (N, 1)
        assert rh_label_arr.ndim == 2 and rh_label_arr.shape[1] == 1, \
            f"rh_label_arr must have shape (N, 1), but got {rh_label_arr.shape}"

        # Convert (N, 1) arrays to (N,)
        lh_label_arr = lh_label_arr.ravel()
        rh_label_arr = rh_label_arr.ravel()

        # Save left & right hemisphere labels
        np.save(lh_label_npy, lh_label_arr)
        np.save(rh_label_npy, rh_label_arr)
    else:
        print(f"WARNING: {label_mat_file} does NOT exist!")


def stack_lh_rh_label(lh_arr, rh_arr):
    """Vertically stack left- and right-hemisphere label arrays.

    Args:
        lh_arr (np.ndarray): Label array for the left hemisphere.
        rh_arr (np.ndarray): Label array for the right hemisphere.

    Returns:
        np.ndarray: Concatenated array containing both hemispheres' labels.
    """
    return np.concatenate((lh_arr, rh_arr))


def compare_Kong2025_MSHBM_apply_medial_wall_mask(label_1_npy, label_2_npy):
    """Remove medial wall vertices (label value 0) from a pair of label arrays.

    Args:
        label_1_npy (np.ndarray): First (e.g. ground-truth) label array.
        label_2_npy (np.ndarray): Second (e.g. simulated) label array.

    Returns:
        tuple: A tuple `(label_1_npy_filtered, label_2_npy_filtered)` of
        np.ndarray, containing `label_1_npy` and `label_2_npy` with any
        vertex that is medial wall (label 0) in either array removed.
    """
    # Create a mask for the medial wall (value 0)
    medialwall_mask = (label_1_npy == 0) | (label_2_npy == 0)

    # Apply the mask to both arrays and remove the medial wall
    label_1_npy_filtered = label_1_npy[~medialwall_mask]
    label_2_npy_filtered = label_2_npy[~medialwall_mask]

    return label_1_npy_filtered, label_2_npy_filtered


def compare_Kong2025_MSHBM_apply_salience_dorsal_attention_mask(
        label_1_npy, label_2_npy):
    """Select Salience and Dorsal Attention Network vertices from label arrays.

    Filters both input label arrays down to vertices whose label belongs to
    the Salience (A&B) and Dorsal Attention (A&B) networks (label values 5,
    15, 9, and 1), keeping only vertices where both arrays fall in this set.

    Args:
        label_1_npy (np.ndarray): First (e.g. ground-truth) label array.
        label_2_npy (np.ndarray): Second (e.g. simulated) label array.

    Returns:
        tuple: A tuple `(label_1_npy_filtered, label_2_npy_filtered)` of
        np.ndarray, containing `label_1_npy` and `label_2_npy` filtered down
        to the Salience/Dorsal Attention network vertices present in both.
    """
    # Selected labels for Salience and Dorsal Attention Networks (A&B)
    mask_label = [5, 15, 9, 1]

    # Create masks for the selected labels in both label_1_npy and label_2_npy
    label_1_mask = np.isin(label_1_npy, mask_label)
    label_2_mask = np.isin(label_2_npy, mask_label)

    # Combine the masks to ensure both labels are in the selected set
    combined_mask = label_1_mask & label_2_mask

    # Apply the combined mask to filter both arrays
    label_1_npy_filtered = label_1_npy[combined_mask]
    label_2_npy_filtered = label_2_npy[combined_mask]

    return label_1_npy_filtered, label_2_npy_filtered


def compare_Kong2025_MSHBM_compute_dice(labels1, labels2, n_labels=17):
    """Compute the mean Dice coefficient between two parcellation label arrays.

    For each label value from 1 to `n_labels`, computes the Dice coefficient
    between the corresponding boolean masks of `labels1` and `labels2`
    (NaN if the label is absent from both arrays), then averages across
    labels while ignoring NaNs.

    Args:
        labels1 (np.ndarray): First parcellation label array.
        labels2 (np.ndarray): Second parcellation label array.
        n_labels (int, optional): Number of parcel labels to score, assumed
            to be numbered 1 to `n_labels`. Defaults to 17.

    Returns:
        float: Mean Dice coefficient across all labels (labels with no
        vertices in either array are excluded from the average).
    """

    dice_scores = []
    for k in range(1, n_labels + 1):
        # Boolean masks of vertices assigned to parcel k in each array
        mask1 = labels1 == k
        mask2 = labels2 == k
        intersection = np.sum(mask1 & mask2)
        denom = np.sum(mask1) + np.sum(mask2)
        # Dice is undefined (NaN) when parcel k is absent from both arrays
        dice_k = (2 * intersection / denom) if denom > 0 else np.nan
        dice_scores.append(dice_k)
    return np.nanmean(dice_scores)


def compare_Kong2025_MSHBM_wrapper(main_config_file, pl_config_file, list_file,
                                   network_type, mode):
    """Compute MSHBM parcellation Dice coefficients for all subjects and sim runs.

    For each subject in `list_file`, converts the ground-truth MSHBM
    parcellation .mat file to .npy and loads/stacks the ground-truth lh/rh
    labels. Then, for every simulation run listed in that subject's
    `sim_run_list_3_{mode}.txt` (deduplicated to one run per simulation
    session and FD/DV threshold combination), converts the simulated MSHBM
    parcellation .mat file to .npy, loads/stacks the simulated lh/rh labels
    (if present), applies the mask appropriate to `network_type`, and saves
    the resulting Dice coefficient as `dice.npy`.

    Args:
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV listing the selected dataset/subID
            combinations to process.
        network_type (str): Parcellation granularity to compare. Must be
            '17_network' (medial wall masked, Dice computed and saved) or
            '4_network' (Salience/Dorsal Attention masked; this path is not
            yet implemented and raises a ValueError).
        mode (str): Simulation selection mode (e.g. 'grid' or 'slide') used
            to locate the `sim_run_list_3_{mode}.txt` file.
    """
    # Read main config
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    pl_name = str(pl_config['pipeline']['pl_name'])
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    MSHBM_w = pl_config['MSHBM']['MSHBM_w']
    MSHBM_MRF = pl_config['MSHBM']['MSHBM_MRF']

    ##########################################################################
    # Ground truth mat2npy
    ##########################################################################
    # Init output_type
    output_type = "MSHBM"

    # Loop through selected subID list
    selected_subID_df = pd.read_csv(list_file,
                                    dtype={
                                        'rm_run_th': str,
                                        'bld': str,
                                        'sim_bld': str,
                                        'highest_sim_bld': str,
                                        'split_bld': str,
                                        'sim_bld': str
                                    })
    for idx, row in selected_subID_df.iterrows():
        # Init run_type
        run_type = "gt"

        # Define identifiers
        dataset, subID = row[['dataset', 'subID']]

        # Define input dirs
        TMS_Kong2025_dir = util.get_TMS_dir_path(run_type, main_config_file,
                                                 dataset, pl_name, subID)
        label_mat_dir = os.path.join(TMS_Kong2025_dir, "MSHBM",
                                     "ind_parcellation", "test_set")

        # Define output dirs
        label_npy_dir = util.get_compare_Kong2025_output_dir_path(
            main_config_file=main_config_file,
            pl_config_file=pl_config_file,
            output_type=output_type,
            dataset=dataset,
            subID=subID,
            sub_output_type=network_type)

        # Convert .mat to .npy files
        compare_Kong2025_MSHBM_mat2npy(pl_config_file, label_mat_dir,
                                       label_npy_dir)

        # Load ground truth MSHBM labels
        lh_label_gt_arr = np.load(
            os.path.join(label_npy_dir,
                         f"lh_labels_w{MSHBM_w}_MRF{MSHBM_MRF}.npy"))
        rh_label_gt_arr = np.load(
            os.path.join(label_npy_dir,
                         f"rh_labels_w{MSHBM_w}_MRF{MSHBM_MRF}.npy"))
        label_gt_arr = stack_lh_rh_label(lh_label_gt_arr, rh_label_gt_arr)

        ##########################################################################
        # Simulation mat2npy
        ##########################################################################
        # Init run_type
        run_type = "ngt"

        # Loop through sim_run-level list
        sim_run_list_file = os.path.join(rep_sel_subID_sim_dir, "result",
                                         "summary", f"sim_run_list_3_{mode}",
                                         dataset, subID,
                                         f"sim_run_list_3_{mode}.txt")
        sim_run_df = pd.read_csv(sim_run_list_file,
                                 dtype={
                                     'rm_run_th': str,
                                     'bld': str,
                                     'sim_bld': str,
                                     'highest_sim_bld': str,
                                     'split_bld': str
                                 })

        # Within each simulation session, only keep 1
        # simulation run (per FD & DV threshold).
        sim_run_df = sim_run_df.drop_duplicates(subset=[
            'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
            f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
        ])

        # Loop through each simulation session
        for _, sim_run_row in sim_run_df.iterrows():
            # Define identifiers
            dataset, subID, rm_run_th, combo, sim_sesID, vary_FD_th, vary_DV_th = sim_run_row[
                [
                    'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
                    f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
                ]]

            # Define input dirs
            TMS_Kong2025_dir = util.get_TMS_dir_path(run_type,
                                                     main_config_file,
                                                     dataset,
                                                     pl_name,
                                                     subID,
                                                     rm_run_th=rm_run_th,
                                                     combo=combo,
                                                     sim_sesID=sim_sesID,
                                                     FD_metric=FD_metric,
                                                     vary_FD_th=vary_FD_th,
                                                     DV_metric=DV_metric,
                                                     vary_DV_th=vary_DV_th)
            label_mat_dir = os.path.join(TMS_Kong2025_dir, "MSHBM",
                                         "ind_parcellation", "test_set")

            # Define output dirs
            label_npy_dir = util.get_compare_Kong2025_output_dir_path(
                main_config_file=main_config_file,
                pl_config_file=pl_config_file,
                output_type=output_type,
                dataset=dataset,
                subID=subID,
                sub_output_type=network_type,
                rm_run_th=rm_run_th,
                combo=combo,
                sim_sesID=sim_sesID,
                FD_metric=FD_metric,
                vary_FD_th=vary_FD_th,
                DV_metric=DV_metric,
                vary_DV_th=vary_DV_th)

            # Convert .mat to .npy files
            compare_Kong2025_MSHBM_mat2npy(pl_config_file, label_mat_dir,
                                           label_npy_dir)

            # Compute and save Dice
            lh_label_ngt_npy = os.path.join(
                label_npy_dir, f"lh_labels_w{MSHBM_w}_MRF{MSHBM_MRF}.npy")
            rh_label_ngt_npy = os.path.join(
                label_npy_dir, f"rh_labels_w{MSHBM_w}_MRF{MSHBM_MRF}.npy")
            if os.path.exists(lh_label_ngt_npy) and os.path.exists(
                    rh_label_ngt_npy):
                # Load simulation labels
                lh_label_ngt_arr = np.load(
                    os.path.join(label_npy_dir,
                                 f"lh_labels_w{MSHBM_w}_MRF{MSHBM_MRF}.npy"))
                rh_label_ngt_arr = np.load(
                    os.path.join(label_npy_dir,
                                 f"rh_labels_w{MSHBM_w}_MRF{MSHBM_MRF}.npy"))
                label_ngt_arr = stack_lh_rh_label(lh_label_ngt_arr,
                                                  rh_label_ngt_arr)

                # Apply masks
                if network_type == "17_network":
                    label_gt_mask_arr, label_ngt_mask_arr = compare_Kong2025_MSHBM_apply_medial_wall_mask(
                        label_gt_arr, label_ngt_arr)

                    # Compute Dice
                    dice = compare_Kong2025_MSHBM_compute_dice(
                        label_gt_mask_arr, label_ngt_mask_arr)
                    np.save(os.path.join(label_npy_dir, "dice.npy"), dice)
                elif network_type == "4_network":
                    label_gt_mask_arr, label_ngt_mask_arr = compare_Kong2025_MSHBM_apply_salience_dorsal_attention_mask(
                        label_gt_arr, label_ngt_arr)
                    raise ValueError(
                        f"`compare_Kong2025_MSHBM_compute_dice` not ready for `4_network` setting"
                    )
                else:
                    raise ValueError(
                        f"Invalid network_type '{network_type}'. Accepted values are '17_network' or '4_network'."
                    )


def submit_compare_Kong2025_MSHBM_wrapper_job(timestamp, jobprefix,
                                              main_config_file, pl_config_file,
                                              list_file, network_type, mode):
    """Submit a PBS job that runs `compare_Kong2025_MSHBM_wrapper`.

    Args:
        timestamp (str): Timestamp string used to namespace the job's
            stdout/stderr directory.
        jobprefix (str): Prefix used to build the PBS job name.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV listing the selected dataset/subID
            combinations to process.
        network_type (str): Parcellation granularity to compare ('17_network'
            or '4_network').
        mode (str): Simulation selection mode (e.g. 'grid' or 'slide').
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rep_compare_Kong2025_dir = pl_config['path'][
        'rep_compare_Kong2025_dep_dir']

    # Determine job std dir & flag file
    std_dir = os.path.join(rep_compare_Kong2025_dir, "MSHBM_wrapper",
                           timestamp, "job")
    os.makedirs(std_dir, exist_ok=True)
    flag_file = os.path.join(std_dir, "done.flag")

    # Construct command
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; sys.path.append(\\\"{TAN2026_TMSMOTION_DIR}/compare\\\"); "
        f"import CBIG_TMSmotion_compare_Kong2025 as compare_Kong2025; "
        f"compare_Kong2025.compare_Kong2025_MSHBM_wrapper(\\\"{main_config_file}\\\", \
            \\\"{pl_config_file}\\\", \\\"{list_file}\\\", \\\"{network_type}\\\", \
            \\\"{mode}\\\")'; "
        f"touch {flag_file}; "
        f"conda deactivate")

    # Init PBS job parameters
    jobname = f"{jobprefix}"
    stderr_file = os.path.join(std_dir, f"{jobprefix}.err")
    stdout_file = os.path.join(std_dir, f"{jobprefix}.out")

    # Submit job to headnode
    util.submit_job(walltime="24:00:00",
                    memory="4G",
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def check_compare_Kong2025_MSHBM_wrapper_job(timestamp, try_count,
                                             main_config_file, pl_config_file):
    """Check whether the `compare_Kong2025_MSHBM_wrapper` job completed.

    Args:
        timestamp (str): Timestamp string identifying the job run, used to
            locate its stdout/stderr directory.
        try_count (int): Retry attempt count. Not used directly in this
            function; kept for interface consistency with other
            `check_*_job` functions.
        main_config_file (str): Path to the main pipeline config JSON file.
        pl_config_file (str): Path to the pipeline-specific config JSON file.

    Returns:
        bool: True if the job's flag file indicates successful completion,
        False otherwise.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rep_compare_Kong2025_dir = pl_config['path'][
        'rep_compare_Kong2025_dep_dir']

    # Check flag file
    flag_file = os.path.join(rep_compare_Kong2025_dir, "MSHBM_wrapper",
                             timestamp, "job", "done.flag")
    return util.check_job_output(flag_file)


#############################################################################
# ===========================================================================
# SUMMARIZE COMPARISON OF OUTPUT TYPES IN SINGLE DATAFRAME
# ===========================================================================
#############################################################################
def filter_compare_num_sim_run_per_sim_ses(pl_config_file, input_df,
                                           output_dir):
    """Keep only a fixed number of simulation sessions per subject/threshold group.

    Within each group of dataset, subID, rm_run_th, and combo, retains the
    first `compare_num_sim_run_per_sim_ses` (from the pipeline config)
    `sim_sesID` values that are also present across every combination of
    varied FD and DV thresholds in that group. The filtered rows are saved
    to `compare_summary_list_3.txt` in `output_dir`.

    Args:
        pl_config_file (str): Path to the pipeline-specific config JSON
            file, used to look up `FD_metric`, `DV_metric`, and
            `compare_num_sim_run_per_sim_ses`.
        input_df (pd.DataFrame): DataFrame of simulation-run rows to filter.
        output_dir (str): Directory to save `compare_summary_list_3.txt`.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    compare_num_sim_run_per_sim_ses = pl_config['compare'][
        'compare_num_sim_run_per_sim_ses']

    # Define the output file path
    output_csv_path = os.path.join(output_dir, "compare_summary_list_3.txt")

    # Columns used for grouping
    group_cols = ['dataset', 'subID', 'rm_run_th', 'combo']
    vary_cols = [f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"]

    # Create an empty list to store the filtered data
    filtered_rows = []

    # Group by dataset, subID, rm_run_th, combo
    for group_keys, group_df in input_df.groupby(group_cols):

        # Get all unique sim_sesID values in the group
        unique_sim_sesIDs = group_df['sim_sesID'].unique()

        # Find sim_sesID values that are present across all vary_{FD_metric}_th and vary_{DV_metric}_th
        valid_sim_sesIDs = []
        for sim_sesID in unique_sim_sesIDs:
            subset = group_df[group_df['sim_sesID'] == sim_sesID]

            # Get the unique pairs of vary_{FD_metric}_th and vary_{DV_metric}_th
            unique_vary_combos = subset[vary_cols].drop_duplicates()

            # Check if all combinations are present
            if len(unique_vary_combos) == len(
                    group_df[vary_cols].drop_duplicates()):
                valid_sim_sesIDs.append(sim_sesID)

        # Select the first compare_num_sim_run_per_sim_ses sim_sesIDs
        selected_sim_sesIDs = valid_sim_sesIDs[:
                                               compare_num_sim_run_per_sim_ses]

        # Filter the DataFrame based on selected sim_sesIDs
        filtered_group = group_df[group_df['sim_sesID'].isin(
            selected_sim_sesIDs)]

        # Append the filtered rows to the list
        filtered_rows.append(filtered_group)

    # Concatenate all filtered groups into a final DataFrame
    output_df = pd.concat(filtered_rows, ignore_index=True)

    # Save the final DataFrame
    output_df.to_csv(output_csv_path, index=False)

    print(f"Filtered output saved to {output_csv_path}")


def generate_compare_summary_subID_list_file(pl_config_file, mode, dataset,
                                             subID, output_dir, TMS_type):
    """Build a per-subject comparison summary combining TMS and MSHBM results.

    Loads the subject's `sim_run_list_3_{mode}.txt`, adds derived columns for
    per-run percent censored frames, FD/DV motion metrics, and whether
    simulation run 1/run 2 are present (validated with assertions), then
    deduplicates to one row per simulation session and FD/DV threshold
    combination. For each remaining row, loads the corresponding
    `eucli_dist.npy` (TMS target distance) and `dice.npy` (17-network MSHBM
    Dice) outputs and appends them as columns (`None` if the file is
    missing). Saves the full result as `compare_summary_list_1.txt`, the
    subset with no missing TMS/MSHBM values as `compare_summary_list_2.txt`,
    and finally calls `filter_compare_num_sim_run_per_sim_ses` to produce
    `compare_summary_list_3.txt`.

    Args:
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        mode (str): Simulation selection mode (e.g. 'grid' or 'slide') used
            to locate the `sim_run_list_3_{mode}.txt` file.
        dataset (str): Dataset identifier for the subject being summarized.
        subID (str): Subject identifier being summarized.
        output_dir (str): Directory to save the summary CSV files.
        TMS_type (str): TMS type identifier used to locate the TMS/MSHBM
            output directories. Currently only 'Kong2025_dep' is a valid
            value.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    num_sim_run_per_sim_ses = pl_config['select_subID']['sim'][
        'num_sim_run_per_sim_ses']

    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']

        # Init output type dict with column names
        output_type_dict = {
            os.path.join(rep_compare_Kong2025_dir, "TMS"):
            ("eucli_dist.npy", "eucli_dist"),
            os.path.join(rep_compare_Kong2025_dir, "MSHBM", "17_network"):
            ("dice.npy", "17_network_dice")
        }
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )

    # Collect rows in a list
    output_rows = []

    # Loop through each simulation run for each threshold
    sim_run_list_file = os.path.join(rep_sel_subID_sim_dir, "result",
                                     "summary", f"sim_run_list_3_{mode}",
                                     dataset, subID,
                                     f"sim_run_list_3_{mode}.txt")
    sim_run_df = pd.read_csv(sim_run_list_file,
                             dtype={
                                 'rm_run_th': str,
                                 'bld': str,
                                 'sim_bld': str,
                                 'highest_sim_bld': str,
                                 'split_bld': str,
                                 'sim_bld': str
                             })

    # Add columns for % censored frames for
    # simulation run 1 & simulation run 2
    pivoted = pd.pivot_table(
        sim_run_df,
        values='sim_run_percent_censored',
        index=['dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID'],
        columns='sim_runID',
        aggfunc='first').reset_index()
    pivoted = pivoted.rename(
        columns={
            'run-01': 'sim_run_percent_censored_run_1',
            'run-02': 'sim_run_percent_censored_run_2'
        })
    sim_run_df = sim_run_df.merge(
        pivoted,
        how='left',
        on=['dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID'])
    # Check that if run-XX present, then must have % censored for that run
    assert sim_run_df.loc[sim_run_df['sim_runID'] == 'run-01', 'sim_run_percent_censored_run_1'].notna().all(), \
        "ERROR: Some run-01 rows have NaN in `sim_run_percent_censored_run_1`"
    if num_sim_run_per_sim_ses > 1:
        assert sim_run_df.loc[sim_run_df['sim_runID'] == 'run-02', 'sim_run_percent_censored_run_2'].notna().all(), \
            "ERROR: Some run-02 rows have NaN in `sim_run_percent_censored_run_2`"
    # Check no NaNs and only allowed values
    allowed = ['run-01', 'run-02']
    assert sim_run_df['sim_runID'].notna().all(
    ), "ERROR: `sim_runID` contains NaN values"
    assert sim_run_df['sim_runID'].isin(allowed).all(
    ), "ERROR: sim_runID contains values other than 'run-01' or 'run-02'"

    # Add columns for FD_metric & DV_metric for
    # simulation run 1 & simulation run 2
    for metric in [FD_metric, DV_metric]:
        pivoted = pd.pivot_table(
            sim_run_df,
            values=f'sim_run_{metric}',
            index=['dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID'],
            columns='sim_runID',
            aggfunc='first').reset_index()
        pivoted = pivoted.rename(
            columns={
                'run-01': f'sim_run_{metric}_run_1',
                'run-02': f'sim_run_{metric}_run_2'
            })
        sim_run_df = sim_run_df.merge(
            pivoted,
            how='left',
            on=['dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID'])
        # Check that if run-XX present, then must have FD/DV_metric for that run
        assert sim_run_df.loc[sim_run_df['sim_runID'] == 'run-01', f'sim_run_{metric}_run_1'].notna().all(), \
            f"ERROR: Some run-01 rows have NaN in `sim_run_{metric}_run_1`"
        if num_sim_run_per_sim_ses > 1:
            assert sim_run_df.loc[sim_run_df['sim_runID'] == 'run-02', f'sim_run_{metric}_run_2'].notna().all(), \
                f"ERROR: Some run-02 rows have NaN in `sim_run_{metric}_run_2`"
        # Check no NaNs and only allowed values
        allowed = ['run-01', 'run-02']
        assert sim_run_df['sim_runID'].notna().all(
        ), "ERROR: `sim_runID` contains NaN values"
        assert sim_run_df['sim_runID'].isin(allowed).all(
        ), "ERROR: sim_runID contains values other than 'run-01' or 'run-02'"

    # Add columns for if simulation run 1 used or not,
    # and if simulation run 2 used or not
    presence = (sim_run_df.assign(flag=1).pivot_table(
        index=['dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID'],
        columns='sim_runID',
        values='flag',
        aggfunc='max',
        fill_value=0).reset_index().rename(columns={
            'run-01': 'sim_run_1',
            'run-02': 'sim_run_2'
        }))
    sim_run_df = sim_run_df.merge(
        presence,
        on=['dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID'],
        how='left')

    # Check that if sim_run_percent_censored_run_X present, then sim_run_X must have 1
    assert sim_run_df.loc[sim_run_df['sim_run_percent_censored_run_1'].notna(), 'sim_run_1'].eq(1).all(), \
        "ERROR: Some rows have sim_run_percent_censored_run_1 not NaN but `sim_run_1` != 1"
    if num_sim_run_per_sim_ses > 1:
        assert sim_run_df.loc[sim_run_df['sim_run_percent_censored_run_2'].notna(), 'sim_run_2'].eq(1).all(), \
            "ERROR: Some rows have sim_run_percent_censored_run_2 not NaN but `sim_run_2` != 1"
    # Check that if sim_run_percent_censored_run_X not present, then sim_run_X must have 0
    assert sim_run_df.loc[sim_run_df['sim_run_percent_censored_run_1'].isna(), 'sim_run_1'].eq(0).all(), \
        "ERROR: Some rows have sim_run_percent_censored_run_1 as NaN but `sim_run_1` != 0"
    if num_sim_run_per_sim_ses > 1:
        assert sim_run_df.loc[sim_run_df['sim_run_percent_censored_run_2'].isna(), 'sim_run_2'].eq(0).all(), \
            "ERROR: Some rows have sim_run_percent_censored_run_2 as NaN but `sim_run_2` != 0"

    # Within each simulation session, only keep 1
    # simulation run (per FD & DV threshold).
    sim_run_df = sim_run_df.drop_duplicates(subset=[
        'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
        f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
    ])

    for _, sim_run_row in sim_run_df.iterrows():
        # Define identifiers
        dataset, subID, rm_run_th, combo, sim_sesID, vary_FD_th, vary_DV_th = sim_run_row[
            [
                'dataset', 'subID', 'rm_run_th', 'combo', 'sim_sesID',
                f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
            ]]

        # Copy row to append new columns
        updated_row = sim_run_row.copy()

        # Loop through each output type
        for base_path, (filename, column_name) in output_type_dict.items():
            # Load output type npy file
            npy_path = os.path.join(
                base_path, dataset, subID, "ngt", rm_run_th, combo, sim_sesID,
                f"{FD_metric}{vary_FD_th}_{DV_metric}{vary_DV_th}", filename)
            try:
                arr = np.load(npy_path)
                # Add array data as a new column
                updated_row[column_name] = arr
            except FileNotFoundError:
                # Handle missing files (if needed)
                updated_row[column_name] = None

        # Collect the updated row
        output_rows.append(updated_row)

    # Create the output DataFrame at once
    output_df_1 = pd.DataFrame(output_rows)

    # Remove rows with missing values in specific columns (the columns from output_type_dict)
    columns_to_check = [
        column_name for _, (filename, column_name) in output_type_dict.items()
    ]
    output_df_2 = output_df_1.dropna(subset=columns_to_check)

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define the path to save the output CSV file
    output_csv_path_1 = os.path.join(output_dir, "compare_summary_list_1.txt")
    output_df_1.to_csv(output_csv_path_1, index=False)

    # Save the cleaned DataFrame (after removing rows with missing values) to a new file
    output_csv_path_2 = os.path.join(output_dir, "compare_summary_list_2.txt")
    output_df_2.to_csv(output_csv_path_2, index=False)

    # Generate `compare_summary_list_3.txt`
    # Keep only first `compare_num_sim_run_per_sim_ses` `sim_sesID`s
    filter_compare_num_sim_run_per_sim_ses(pl_config_file, output_df_2,
                                           output_dir)


def submit_generate_compare_summary_subID_list_file_job(
        timestamp, jobprefix, pl_config_file, mode, dataset, subID, TMS_type):
    """Submit a PBS job that runs `generate_compare_summary_subID_list_file`.

    Args:
        timestamp (str): Timestamp string used to namespace the job's
            stdout/stderr directory.
        jobprefix (str): Prefix used to build the PBS job name (combined
            with `dataset` and `subID`).
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        mode (str): Simulation selection mode (e.g. 'grid' or 'slide').
        dataset (str): Dataset identifier for the subject being summarized.
        subID (str): Subject identifier being summarized.
        TMS_type (str): TMS type identifier, currently only 'Kong2025_dep'
            is a valid value.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )
    # Determine file paths and directories
    output_dir = os.path.join(rep_compare_Kong2025_dir, "summary", dataset,
                              subID)

    # Create job std dir & determine flag file
    std_dir = os.path.join(rep_compare_Kong2025_dir, "summary", timestamp,
                           "job", dataset, subID)
    os.makedirs(std_dir, exist_ok=True)
    flag_file = os.path.join(std_dir, "done.flag")

    # Construct command
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; sys.path.append(\\\"{TAN2026_TMSMOTION_DIR}/compare\\\"); "
        f"import CBIG_TMSmotion_compare_Kong2025 as compare_Kong2025; "
        f"compare_Kong2025.generate_compare_summary_subID_list_file(\\\"{pl_config_file}\\\", \
            \\\"{mode}\\\", \\\"{dataset}\\\", \\\"{subID}\\\", \\\"{output_dir}\\\", \
            \\\"{TMS_type}\\\")'; "
        f"touch {flag_file}; "
        f"conda deactivate")

    # Init PBS job parameters
    jobname = f"{jobprefix}_{dataset}_{subID}"
    stderr_file = os.path.join(std_dir, f"{dataset}_{subID}.err")
    stdout_file = os.path.join(std_dir, f"{dataset}_{subID}.out")

    # Submit job to headnode
    util.submit_job(walltime="24:00:00",
                    memory="2G",
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def submit_all_generate_compare_summary_subID_list_file_job(
        timestamp, jobprefix, pl_config_file, list_file, mode, TMS_type):
    """Submit `generate_compare_summary_subID_list_file` jobs for all subjects.

    Reads the selected dataset/subID combinations from `list_file` and
    submits one `submit_generate_compare_summary_subID_list_file_job` PBS
    job per subject.

    Args:
        timestamp (str): Timestamp string used to namespace each job's
            stdout/stderr directory.
        jobprefix (str): Prefix used to build each job's PBS job name.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV listing the selected dataset/subID
            combinations to process.
        mode (str): Simulation selection mode (e.g. 'grid' or 'slide').
        TMS_type (str): TMS type identifier, currently only 'Kong2025_dep'
            is a valid value.
    """

    # Read list all selected subIDs
    selected_subID_df = pd.read_csv(list_file,
                                    dtype={
                                        'rm_run_th': str,
                                        'bld': str,
                                        'sim_bld': str,
                                        'highest_sim_bld': str,
                                        'split_bld': str,
                                        'highest_bld': str
                                    })

    # Loop through each selected subID
    for _, row in selected_subID_df.iterrows():
        # Init identifiers
        dataset, subID = row[['dataset', 'subID']]

        # Submit job for 1 selected subID
        submit_generate_compare_summary_subID_list_file_job(
            timestamp, jobprefix, pl_config_file, mode, dataset, subID,
            TMS_type)


def check_generate_compare_summary_subID_list_file_job(timestamp, try_count,
                                                       pl_config_file, dataset,
                                                       subID, TMS_type,
                                                       all_job_success,
                                                       fail_list):
    """Check the `generate_compare_summary_subID_list_file` job for one subject.

    Checks the job's flag file to verify successful completion. On failure,
    marks `all_job_success` as False and appends `dataset,subID` to the
    shared `fail_list` CSV (creating it with a header if it does not yet
    exist).

    Args:
        timestamp (str): Timestamp string identifying the job run, used to
            locate its flag file.
        try_count (int): Retry attempt count. Not used directly in this
            function; kept for interface consistency with other
            `check_*_job` functions.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        dataset (str): Dataset identifier for the subject being checked.
        subID (str): Subject identifier being checked.
        TMS_type (str): TMS type identifier, currently only 'Kong2025_dep'
            is a valid value.
        all_job_success (bool): Running success flag accumulated across
            subjects; set to False if this subject's job failed.
        fail_list (str): Path to the shared CSV file that records failed
            dataset/subID pairs for this try.

    Returns:
        tuple: `(all_job_success, fail_list)` — the (possibly updated)
        running success flag and the path to the fail-list CSV.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )

    # Check flag file for success
    flag_file = os.path.join(rep_compare_Kong2025_dir, "summary", timestamp,
                             "job", dataset, subID, "done.flag")
    id_success = util.check_job_output(flag_file)
    if not id_success:
        all_job_success = False
        # Handle failed ID list (shared across all subIDs for this try)
        if not os.path.exists(fail_list):
            with open(fail_list, 'w') as f:
                f.write("dataset,subID\n")
        with open(fail_list, 'a') as f:
            f.write(f"{dataset},{subID}\n")

    return all_job_success, fail_list


def check_all_generate_compare_summary_subID_list_file_job(
        timestamp, try_count, pl_config_file, list_file, TMS_type):
    """Check `generate_compare_summary_subID_list_file` jobs for all subjects.

    Reads the selected dataset/subID combinations from `list_file` and calls
    `check_generate_compare_summary_subID_list_file_job` for each, tracking
    an overall success flag and a shared fail-list CSV (named using
    `try_count`) for this attempt.

    Args:
        timestamp (str): Timestamp string identifying the job run.
        try_count (int): Retry attempt number, used to name the fail-list
            CSV file (`fail_list_try_{try_count}.txt`).
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV listing the selected dataset/subID
            combinations to check.
        TMS_type (str): TMS type identifier, currently only 'Kong2025_dep'
            is a valid value.

    Returns:
        tuple: `(all_job_success, fail_list)` where `all_job_success` (bool)
        is True only if every subject's job succeeded, and `fail_list` (str)
        is the path to the CSV recording any failed dataset/subID pairs.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )

    # Define 1 shared fail list for this try, aggregated across all subIDs
    timestamp_dir = os.path.join(rep_compare_Kong2025_dir, "summary",
                                 timestamp)
    os.makedirs(timestamp_dir, exist_ok=True)
    fail_list = os.path.join(timestamp_dir, f"fail_list_try_{try_count}.txt")

    # Read list all selected subIDs
    selected_subID_df = pd.read_csv(list_file,
                                    dtype={
                                        'rm_run_th': str,
                                        'bld': str,
                                        'sim_bld': str,
                                        'highest_sim_bld': str,
                                        'split_bld': str,
                                        'highest_bld': str
                                    })

    # Init check job vars
    all_job_success = True

    # Loop through each selected subID
    for _, row in selected_subID_df.iterrows():
        # Init identifiers
        dataset, subID = row[['dataset', 'subID']]

        # Check job for 1 selected subID
        all_job_success, fail_list = check_generate_compare_summary_subID_list_file_job(
            timestamp, try_count, pl_config_file, dataset, subID, TMS_type,
            all_job_success, fail_list)

    return all_job_success, fail_list


def concat_all_compare_summary_subID_list_file(pl_config_file,
                                               list_file,
                                               TMS_type='Kong2025_dep'):
    """Concatenate per-subject comparison summaries into a single file.

    Reads the selected dataset/subID combinations from `list_file` and, for
    each, loads that subject's `compare_summary_list_3.txt` (printing a
    warning and skipping it if missing). Concatenates all loaded summaries
    and saves the combined result as `compare_summary_list_4.txt` in the
    comparison output directory.

    Args:
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV listing the selected dataset/subID
            combinations to concatenate.
        TMS_type (str, optional): TMS type identifier. Defaults to
            'Kong2025_dep', currently the only accepted value.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )

    # Define output dirs
    output_dir = os.path.join(rep_compare_Kong2025_dir, "summary")

    # Read list all selected subIDs
    selected_subID_df = pd.read_csv(list_file,
                                    dtype={
                                        'rm_run_th': str,
                                        'bld': str,
                                        'sim_bld': str,
                                        'highest_sim_bld': str,
                                        'split_bld': str,
                                        'sim_bld': str
                                    })

    # Initialize list to collect DataFrames
    all_subID_data = []

    # Loop through each selected subID
    for _, row in selected_subID_df.iterrows():
        # Init identifiers
        dataset, subID = row[['dataset', 'subID']]

        # Define input files
        subID_summary_list = os.path.join(rep_compare_Kong2025_dir, "summary",
                                          dataset, subID,
                                          "compare_summary_list_3.txt")

        try:
            # Try reading the subID summary list CSV
            subID_df = pd.read_csv(subID_summary_list,
                                   dtype={
                                       'rm_run_th': str,
                                       'bld': str,
                                       'sim_bld': str,
                                       'highest_sim_bld': str,
                                       'split_bld': str,
                                       'sim_bld': str
                                   })
            all_subID_data.append(subID_df)  # Add this DataFrame to the list
        except FileNotFoundError:
            # Handle missing file with a warning message
            print(
                f"Warning: The file for subID {subID} in dataset {dataset} does not exist: {subID_summary_list}"
            )

    if all_subID_data:
        # Concatenate all DataFrames (vertically stack them)
        concatenated_df = pd.concat(all_subID_data, ignore_index=True)

        # Define the output path for the concatenated result
        output_csv_path = os.path.join(output_dir,
                                       f"compare_summary_list_4.txt")

        # Save the concatenated DataFrame to CSV
        concatenated_df.to_csv(output_csv_path, index=False)
        print(f"Concatenated data saved to: {output_csv_path}")
    else:
        print(
            "No valid subID summary files were found. No data to concatenate.")


def submit_concat_all_compare_summary_subID_list_file_job(
        timestamp, jobprefix, pl_config_file, list_file, TMS_type):
    """Submit a PBS job that runs `concat_all_compare_summary_subID_list_file`.

    Args:
        timestamp (str): Timestamp string used to namespace the job's
            stdout/stderr directory.
        jobprefix (str): Prefix used to build the PBS job name.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        list_file (str): Path to a CSV listing the selected dataset/subID
            combinations to concatenate.
        TMS_type (str): TMS type identifier, currently only 'Kong2025_dep'
            is a valid value.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )

    # Determine job std dir & flag file
    std_dir = os.path.join(rep_compare_Kong2025_dir, "concat_wrapper",
                           timestamp, "job")
    os.makedirs(std_dir, exist_ok=True)
    flag_file = os.path.join(std_dir, "done.flag")

    # Construct command
    cmd_script = (
        f"source CBIG_init_conda; conda activate CBIG_TMSmotion; "
        f"python -c 'import sys; sys.path.append(\\\"{TAN2026_TMSMOTION_DIR}/compare\\\"); "
        f"import CBIG_TMSmotion_compare_Kong2025 as compare_Kong2025; "
        f"compare_Kong2025.concat_all_compare_summary_subID_list_file(\\\"{pl_config_file}\\\", \
            \\\"{list_file}\\\", TMS_type=\\\"{TMS_type}\\\")'; "
        f"touch {flag_file}; "
        f"conda deactivate")

    # Init PBS job parameters
    jobname = f"{jobprefix}"
    stderr_file = os.path.join(std_dir, f"{jobprefix}.err")
    stdout_file = os.path.join(std_dir, f"{jobprefix}.out")

    # Submit job to headnode
    util.submit_job(walltime="24:00:00",
                    memory="4G",
                    ncpus="1",
                    jobname=jobname,
                    stderr_file=stderr_file,
                    stdout_file=stdout_file,
                    cmd_script=cmd_script)


def check_concat_all_compare_summary_subID_list_file_job(
        timestamp, try_count, pl_config_file, TMS_type):
    """Check whether the `concat_all_compare_summary_subID_list_file` job completed.

    Args:
        timestamp (str): Timestamp string identifying the job run, used to
            locate its flag file.
        try_count (int): Retry attempt count. Not used directly in this
            function; kept for interface consistency with other
            `check_*_job` functions.
        pl_config_file (str): Path to the pipeline-specific config JSON file.
        TMS_type (str): TMS type identifier, currently only 'Kong2025_dep'
            is a valid value.

    Returns:
        bool: True if the job's flag file indicates successful completion,
        False otherwise.
    """

    # Read pipeline config
    pl_config = util.read_json_file(pl_config_file)
    if TMS_type == 'Kong2025_dep':
        rep_compare_Kong2025_dir = pl_config['path'][
            'rep_compare_Kong2025_dep_dir']
    else:
        raise ValueError(
            f"Invalid TMS_type '{TMS_type}'. Accepted value is only 'Kong2025_dep'."
        )

    # Check flag file
    flag_file = os.path.join(rep_compare_Kong2025_dir, "concat_wrapper",
                             timestamp, "job", "done.flag")
    return util.check_job_output(flag_file)


#############################################################################
# ===========================================================================
# MAIN
# ===========================================================================
#############################################################################


def main():
    """Run the full Kong2025 TMS/MSHBM comparison for one replication pipeline.

    Reads the main and pipeline (given as the first command-line argument)
    config files, then runs `compare_Kong2025_TMS_wrapper` (TMS target
    Euclidean distance) and `compare_Kong2025_MSHBM_wrapper` (17-network
    MSHBM Dice) for the subjects listed in `selected_subID_list_2.txt`.
    After `CBIG_TMSmotion_compare_Kong2025_wrapper.py`'s
    `submit_all_generate_compare_summary_subID_list_file_job` has been run
    separately to build per-subject summaries, this also concatenates those
    summaries across subjects via
    `concat_all_compare_summary_subID_list_file`.
    """
    # Read main config
    main_config_file = os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                                    'config', 'main_config.json')
    util.read_json_file(main_config_file)

    # Read pipeline config
    pl_config_file = sys.argv[1]
    pl_config = util.read_json_file(pl_config_file)

    # Define slide or grid mode
    mode = 'grid'

    # Define TMS type
    TMS_type = 'Kong2025_dep'
    rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']

    # Define input files
    list_file = os.path.join(rep_sel_subID_sim_dir, "result", "summary",
                             "selected_subID_list_2.txt")

    ####################################
    # Generate variables of interest
    ####################################

    # run compare TMS wrapper
    compare_Kong2025_TMS_wrapper(main_config_file,
                                 pl_config_file,
                                 list_file,
                                 mode,
                                 TMS_type=TMS_type)

    network_type = "17_network"
    compare_Kong2025_MSHBM_wrapper(main_config_file, pl_config_file, list_file,
                                   network_type, mode)

    # after that, run `CBIG_TMSmotion_compare_Kong2025_wrapper.py`'s
    # `submit_all_generate_compare_summary_subID_list_file_job`.

    ####################################
    # Concatenate variables of interest
    # across subIDs
    ####################################
    # after running `CBIG_TMSmotion_compare_Kong2025_wrapper.py`'s
    # `submit_all_generate_compare_summary_subID_list_file_job`
    list_file = os.path.join(rep_sel_subID_sim_dir, "result", "summary",
                             "selected_subID_list_2.txt")
    concat_all_compare_summary_subID_list_file(pl_config_file,
                                               list_file,
                                               TMS_type=TMS_type)


if __name__ == "__main__":
    main()
