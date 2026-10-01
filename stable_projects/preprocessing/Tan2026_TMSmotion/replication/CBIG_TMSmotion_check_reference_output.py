'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
"""
Checks a replication run's `compare_summary_list_4.txt` against the
reference output for the corresponding pipeline & full/lite mode.

For each unique (dataset, subID, rm_run_th, combo, sim_run_bin, sim_sesID,
vary_FD_th (or vary_FDRMS_th), vary_DVARS_th), checks that `eucli_dist` and
`17_network_dice` are within tolerance of the reference. Also flags any row
present in one output but missing from the other.

Examples:

python $TAN2026_TMSMOTION_DIR/replication/CBIG_TMSmotion_check_reference_output.py pl1 lite

python $TAN2026_TMSMOTION_DIR/replication/CBIG_TMSmotion_check_reference_output.py pl2 full

"""
#############################################################################
# ===========================================================================
# IMPORTS
# ===========================================================================
#############################################################################
import os
import sys
import numpy as np
import pandas as pd

TAN2026_TMSMOTION_DIR = os.environ.get('TAN2026_TMSMOTION_DIR')
sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'util'))
import CBIG_TMSmotion_util as util

CBIG_REPDATA_DIR = os.environ.get('CBIG_REPDATA_DIR')

# Only 'Kong2025_dep' is currently in active use across the replication
# pipeline (matches compare/CBIG_TMSmotion_compare_Kong2025.py's `main()`).
TMS_TYPE_DIR = "compare_Kong2025_dep"

# Columns identifying 1 unique simulated result (vary_FD_th/vary_DV_th are
# appended dynamically, since their name depends on FD_metric/DV_metric)
KEY_COLUMNS = [
    'dataset', 'subID', 'rm_run_th', 'combo', 'sim_run_bin', 'sim_sesID'
]

# Metric columns to compare within tolerance
METRIC_COLUMNS = ['eucli_dist', '17_network_dice']

TOLERANCE = 1e-4


#############################################################################
# ===========================================================================
# FUNCTIONS
# ===========================================================================
#############################################################################
def read_compare_summary_list_4(filepath, key_columns):
    """Read and validate a `compare_summary_list_4.txt` file for comparison.

    Loads the CSV at `filepath`, verifies that all `key_columns` plus the
    module-level `METRIC_COLUMNS` are present, then returns only those
    columns with one row per unique combination of `key_columns` values.

    Args:
        filepath (str): Path to the `compare_summary_list_4.txt` file to
            read.
        key_columns (list): Column names that together uniquely identify
            one simulated result; used to check for duplicates and to
            validate required columns.

    Returns:
        pd.DataFrame: The key and metric columns, deduplicated on
            `key_columns`.

    Raises:
        FileNotFoundError: If `filepath` does not exist.
        ValueError: If any of `key_columns` or `METRIC_COLUMNS` is missing
            from the file.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Output file does NOT exist: {filepath}")

    dtype_dict = {col: str for col in KEY_COLUMNS}
    df = pd.read_csv(filepath, dtype=dtype_dict)

    missing_columns = [
        col for col in key_columns + METRIC_COLUMNS if col not in df.columns
    ]
    if missing_columns:
        raise ValueError(
            f"{filepath} is missing expected column(s): {missing_columns}")

    # Keep only the columns needed for comparison, 1 row per unique key
    return df[key_columns + METRIC_COLUMNS].drop_duplicates(subset=key_columns)


def check_reference_output(pl_name, mode):
    """Compare a replication run's summary output against the reference output.

    Reads the pipeline config for `pl_name`, then loads both the user's
    `compare_summary_list_4.txt` output and the corresponding reference
    output for `mode` ("full" or "lite"). Rows are merged on the key
    columns; any row present in only one of the two files is reported as
    an error, and for matched rows, `eucli_dist` and `17_network_dice`
    are checked to be within `TOLERANCE` of the reference. Prints a
    summary of any discrepancies found. If any check fails, the process
    exits with status code 1 (via `sys.exit`).

    Args:
        pl_name (str): Name of the pipeline whose output to check (e.g.
            "pl1").
        mode (str): Either "full" or "lite"; selects which reference
            output to compare against.
    """
    # Read pipeline config
    pl_config_file = os.path.join(TAN2026_TMSMOTION_DIR, 'replication',
                                  pl_name, 'config', 'pipeline.json')
    pl_config = util.read_json_file(pl_config_file)
    FD_metric = pl_config['preproc']['FD_metric']
    DV_metric = pl_config['preproc']['DV_metric']
    key_columns = KEY_COLUMNS + [
        f"vary_{FD_metric}_th", f"vary_{DV_metric}_th"
    ]

    # Define user & reference output files
    user_file = os.path.join(pl_config['path']['rep_compare_Kong2025_dep_dir'],
                             'summary', 'compare_summary_list_4.txt')
    ref_file = os.path.join(CBIG_REPDATA_DIR, 'stable_projects',
                            'preprocessing', 'Tan2026_TMSmotion', 'ref_result',
                            pl_name, mode, TMS_TYPE_DIR, 'summary',
                            'compare_summary_list_4.txt')

    user_df = read_compare_summary_list_4(user_file, key_columns)
    ref_df = read_compare_summary_list_4(ref_file, key_columns)

    # Merge reference & user outputs on the key columns
    merged_df = ref_df.merge(user_df,
                             on=key_columns,
                             how='outer',
                             suffixes=('_ref', '_user'),
                             indicator=True)

    has_error = False

    # Flag rows present in reference but missing from user output, & vice versa
    missing_in_user = merged_df[merged_df['_merge'] == 'left_only']
    extra_in_user = merged_df[merged_df['_merge'] == 'right_only']
    if not missing_in_user.empty:
        has_error = True
        print(f"ERROR: {len(missing_in_user)} row(s) in reference are "
              f"MISSING from user output:")
        print(missing_in_user[key_columns].to_string(index=False))
    if not extra_in_user.empty:
        has_error = True
        print(f"ERROR: {len(extra_in_user)} row(s) in user output are "
              f"EXTRA (not present in reference):")
        print(extra_in_user[key_columns].to_string(index=False))

    # Check metric columns are within tolerance, for matched rows
    matched_df = merged_df[merged_df['_merge'] == 'both']
    for metric in METRIC_COLUMNS:
        ref_col = f"{metric}_ref"
        user_col = f"{metric}_user"
        is_close = np.isclose(matched_df[ref_col],
                              matched_df[user_col],
                              atol=TOLERANCE,
                              rtol=0,
                              equal_nan=True)
        out_of_tol = matched_df[~is_close]
        if not out_of_tol.empty:
            has_error = True
            print(f"ERROR: {len(out_of_tol)} row(s) have '{metric}' "
                  f"differing from reference by more than {TOLERANCE}:")
            report_df = out_of_tol[key_columns + [ref_col, user_col]].copy()
            report_df['abs_diff'] = (report_df[ref_col] -
                                     report_df[user_col]).abs()
            print(report_df.to_string(index=False))

    if has_error:
        print(f"Replication FAILED for {pl_name} ({mode}).")
        sys.exit(1)
    else:
        print(f"Replication success for {pl_name}")


#############################################################################
# ===========================================================================
# MAIN
# ===========================================================================
#############################################################################
def main():
    """Validate command-line arguments and run the reference-output check.

    Reads the pipeline name and mode from the command line
    (`sys.argv[1]` and `sys.argv[2]`), validates that mode is either
    "full" or "lite", then calls `check_reference_output`.

    Raises:
        ValueError: If `mode` is not "full" or "lite".
    """
    pl_name = sys.argv[1]
    mode = sys.argv[2]
    if mode not in ('full', 'lite'):
        raise ValueError(
            f"Invalid mode '{mode}'. Accepted values are 'full' or 'lite'.")

    check_reference_output(pl_name, mode)


if __name__ == "__main__":
    main()
