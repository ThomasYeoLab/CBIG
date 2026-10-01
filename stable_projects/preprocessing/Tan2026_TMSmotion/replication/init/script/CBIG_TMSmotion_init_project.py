'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

Examples:

python $TAN2026_TMSMOTION_DIR/replication/init/script/CBIG_TMSmotion_init_project.py

'''

# ===========================================================================
# IMPORTS
# ===========================================================================
import os
import sys
import getpass

TAN2026_TMSMOTION_DIR = os.environ.get('TAN2026_TMSMOTION_DIR')
sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'util'))
CBIG_REPDATA_DIR = os.environ.get('CBIG_REPDATA_DIR')
TAN2026_TMSMOTION_DATA_DIR = os.environ.get('TAN2026_TMSMOTION_DATA_DIR')
import CBIG_TMSmotion_util as util


# ===========================================================================
# FUNCTIONS
# ===========================================================================
def generate_main_config():
    """Generate the main configuration file shared across all pipelines and subIDs.

    Builds a config dictionary covering the job scheduler username,
    shared paths (e.g. `first_N_T1_list`, `data_dir`, `rep_dir`), and the
    struct-preprocessing `reconall_map`, then writes it to
    `<TAN2026_TMSMOTION_DIR>/replication/config/main_config.json`.

    For slice-timing-correction (STC)-related information, the possible
    meanings are as follows:

    For datasets with JSON provided:
        'one_dataset_one_json' = 1 dataset has 1 JSON
        'one_run_one_json' = 1 run has 1 JSON

    For datasets without JSON provided, documentation suggests:
        'sequential_increasing' = slices obtained inferior -> superior
            sequentially, like 1, 2, ..., N
        'alternating_increasing' = interleaved slices obtained inferior
            -> superior, like 1, 3, 5, ..., 2, 4, 6
        'alternating_increasing_2' = interleaved slices obtained
            inferior -> superior, like 2, 4, 6, ..., 1, 3, 5
    """

    # Define configuration dictionary
    config = {
        "jobscheduler": {
            "username": f"{getpass.getuser()}"
        },
        "path": {
            "first_N_T1_list":
            f"{CBIG_REPDATA_DIR}/stable_projects/preprocessing/Tan2026_TMSmotion/input_data/first_N_T1_list.txt",
            "data_dir": f"{TAN2026_TMSMOTION_DATA_DIR}",
            "rep_dir": f"{TAN2026_TMSMOTION_DIR}/replication"
        },
        "preproc_struct": {
            "reconall_map": [{
                "dataset": "Newbold",
                "subID": "sub-cast1",
                "map2": {
                    "dataset": "MSC",
                    "subID": "sub-MSC02"
                }
            }]
        }
    }

    # Define the output file name
    output_filename = f"{TAN2026_TMSMOTION_DIR}/replication/config/main_config.json"
    util.write_config_dict2json(config, output_filename)


def main():
    """Generate and write the main configuration file.

    Thin wrapper that calls `generate_main_config`.
    """
    # Generate main config
    generate_main_config()


if __name__ == "__main__":
    main()
