'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

This script is called by examples/CBIG_TMSmotion_examples_wrapper.sh (not meant
to be run directly by the user), e.g.:

Examples:

bash $TAN2026_TMSMOTION_DIR/examples/CBIG_TMSmotion_examples_wrapper.sh
'''

#############################################################################
# IMPORTS
#############################################################################
import os
import sys
import numpy as np
import pandas as pd
import warnings

warnings.simplefilter(action="ignore",
                      category=pd.errors.SettingWithCopyWarning)


#############################################################################
# MAIN
#############################################################################
def main():
    '''Simulate rs-fMRI runs from example non-ground-truth data for two participants.

    For each of the two synthetic participants ("sub-01" and "sub-02"),
    loads their windowed motion QC records, simulates candidate runs via
    `CBIG_TMSmotion_select_subID_sim.simulate_run` across a set of
    sample bins, and saves the selected window indices to a `.npy` file
    (shape: number of windows per run x number of simulated runs) named
    `<subID>.npy` in `<TAN2026_TMSMOTION_DIR>/<out_name>/output/`.

    The output subdirectory name (`out_name`, under
    `TAN2026_TMSMOTION_DIR`) is read from `sys.argv[1]`.

    Raises:
        EnvironmentError: If the `TAN2026_TMSMOTION_DIR` environment
            variable is not set.
    '''

    # --------------------------------------------------------------------------
    # Environment setup
    # --------------------------------------------------------------------------
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
    sys.path.append(os.path.join(PROJECT_DIR, "util"))
    import CBIG_TMSmotion_util as util
    sys.path.append(os.path.join(PROJECT_DIR, "examples", "scripts"))
    import CBIG_TMSmotion_select_subID_sim_examples as select_subID_sim

    # ----------------------------------------------------------------------
    # Config
    # ----------------------------------------------------------------------
    pl_config_file = os.path.join(PROJECT_DIR, "examples", "input",
                                  "pipeline.json")

    windowID_list = os.path.join(PROJECT_DIR, "examples", "input",
                                 "window_list_1.txt")
    out_name = sys.argv[1]
    output_dir = os.path.join(PROJECT_DIR, out_name, "output")
    os.makedirs(output_dir, exist_ok=True)

    # ----------------------------------------------------------------------
    # Sampling parameters
    # ----------------------------------------------------------------------
    dataset = "examples_dataset"
    subID_list = ["sub-01", "sub-02"]
    sample_bin_interval = 10
    sample_bin_min = -10
    sample_bin_max = 100
    sample_bin_list = select_subID_sim.generate_bin_list(
        sample_bin_min, sample_bin_max, sample_bin_interval)
    target_count_per_bin = 100

    # ----------------------------------------------------------------------
    # Load dataframe containing window
    # ----------------------------------------------------------------------
    df = pd.read_csv(windowID_list, dtype={"bld": str})

    # ----------------------------------------------------------------------
    # Run simulation for each participant
    # ----------------------------------------------------------------------
    for subID in subID_list:

        print(f"Running simulation for {subID}...")

        df_subID = util.filter_dataset_subID_from_all_subID_df(
            df, dataset, subID)
        df_subID = df_subID.reset_index(drop=True)

        selected_window_idx, _, _, _, _ = (select_subID_sim.simulate_run(
            pl_config_file,
            df_subID,
            sample_bin_list,
            target_count_per_bin,
            dataset,
            subID,
            examples=True))

        # Save numpy array
        output_file = os.path.join(output_dir, f"{subID}.npy")
        np.save(output_file, selected_window_idx)

    print("Examples simulations completed successfully.\n")


if __name__ == "__main__":
    main()
