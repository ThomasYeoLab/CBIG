'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
"""
Examples:

python $TAN2026_TMSMOTION_DIR/replication/TMS_Kong2025/script/CBIG_TMSmotion_Kong2025_gt_wrapper.py \
    $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json

python $TAN2026_TMSMOTION_DIR/replication/TMS_Kong2025/script/CBIG_TMSmotion_Kong2025_gt_wrapper.py \
    $TAN2026_TMSMOTION_DIR/replication/pl2/config/pipeline.json

"""
#############################################################################
# ===========================================================================
# IMPORTS
# ===========================================================================
#############################################################################

import os
import sys
import pandas as pd
from datetime import datetime

TAN2026_TMSMOTION_DIR = os.environ.get('TAN2026_TMSMOTION_DIR')
sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'util'))
import CBIG_TMSmotion_util as util

CBIG_CODE_DIR = util.get_cbig_code_dir()

sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'TMS', 'Kong2025'))
import CBIG_TMSmotion_Kong2025 as TMS_Kong2025

#############################################################################
# ===========================================================================
# GENERAL FOR ALL JOBS
# ===========================================================================
#############################################################################
# Read main config
main_config_file = os.path.join(TAN2026_TMSMOTION_DIR, 'replication', 'config',
                                'main_config.json')
main_config = util.read_json_file(main_config_file)
data_dir = main_config['path']['data_dir']
username = str(main_config['jobscheduler']['username'])

# Read pipeline config
pl_config_file = sys.argv[1]
pl_config = util.read_json_file(pl_config_file)
pl_name = str(pl_config['pipeline']['pl_name'])
rep_sel_subID_gt_dir = pl_config['path']['rep_sel_subID_gt_dir']
rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']
rep_TMS_Kong2025_gt_gt_run_list_2_dir = pl_config['path'][
    'rep_TMS_Kong2025_gt_gt_run_list_2_dir']
# Optional 2nd CLI arg selects the TMS target type; defaults to 'dep' if omitted/empty
target_type = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] else 'dep'

# Init convenience vars
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
run_type = "gt"
reduce_mem_flag = "1"

# Prepare subID-level lists
all_selected_subID_df = pd.read_csv(
    os.path.join(rep_sel_subID_sim_dir, "result", "summary",
                 "selected_subID_list_2.txt"),
    dtype={
        'bld': str,
        'sim_bld': str,
        'highest_split_bld': str
    })
subID_list = list(
    all_selected_subID_df.groupby(['dataset', 'subID'], sort=False))
for (dataset, subID), _ in all_selected_subID_df.groupby(['dataset', 'subID']):
    #############################################################################
    # ===========================================================================
    # GT_LIST_2 GENERATION
    # ===========================================================================
    #############################################################################
    # Init preproc job-specific vars
    check_job_func = TMS_Kong2025.check_generate_gt_list_2_job
    submit_job_func = TMS_Kong2025.submit_generate_gt_list_2_job
    try_count = 1
    try_max = 2  # inclusive of 1st try
    jobprefix = f"GGl_{pl_name}"

    check_kwargs = {
        'timestamp': timestamp,
        'try_count': try_count,
        'main_config_file': main_config_file,
        'pl_config_file': pl_config_file,
        'dataset': dataset,
        'subID': subID,
        'run_type': run_type,
        'reduce_mem_flag': reduce_mem_flag
    }
    submit_kwargs = {
        'timestamp': timestamp,
        'jobprefix': jobprefix,
        'main_config_file': main_config_file,
        'pl_config_file': pl_config_file,
        'dataset': dataset,
        'subID': subID,
        'run_type': run_type,
        'reduce_mem_flag': reduce_mem_flag
    }

    # Check & re-submit job until hit try limit
    util.check_and_resubmit_job(check_job_func,
                                submit_job_func,
                                timestamp,
                                try_count,
                                try_max,
                                jobprefix,
                                username,
                                fail_list_flag=False,
                                check_kwargs=check_kwargs,
                                submit_kwargs=submit_kwargs)

    #############################################################################
    # ===========================================================================
    # SPLIT GROUND TRUTH RUNS INTO SUB-RUNS
    # ===========================================================================
    #############################################################################
    # Init preproc job-specific vars
    gt_list_2 = os.path.join(rep_TMS_Kong2025_gt_gt_run_list_2_dir,
                             "gt_list_2", dataset, subID, "gt_list_2.txt")
    check_job_func = TMS_Kong2025.check_all_gt_split_job
    submit_job_func = TMS_Kong2025.submit_all_gt_split_job
    try_count = 1
    try_max = 2  # inclusive of 1st try
    jobprefix = f"GSp_{pl_name}"
    check_kwargs = {
        'timestamp': timestamp,
        'try_count': try_count,
        'main_config_file': main_config_file,
        'pl_config_file': pl_config_file,
        'list_file': gt_list_2,
        'run_type': run_type,
        'reduce_mem_flag': reduce_mem_flag
    }
    submit_kwargs = {
        'timestamp': timestamp,
        'jobprefix': jobprefix,
        'main_config_file': main_config_file,
        'pl_config_file': pl_config_file,
        'run_type': run_type,
        'list_file': gt_list_2,
        'reduce_mem_flag': reduce_mem_flag
    }

    # Check & re-submit job until hit try limit
    util.check_and_resubmit_job(check_job_func,
                                submit_job_func,
                                timestamp,
                                try_count,
                                try_max,
                                jobprefix,
                                username,
                                fail_list_flag=True,
                                check_kwargs=check_kwargs,
                                submit_kwargs=submit_kwargs)

    #############################################################################
    # ===========================================================================
    # TMS KONG2025
    # ===========================================================================
    #############################################################################
    gt_list_2 = os.path.join(rep_TMS_Kong2025_gt_gt_run_list_2_dir,
                             "gt_list_2", dataset, subID, "gt_list_2.txt")
    # Init preproc job-specific vars
    check_job_func = TMS_Kong2025.check_TMS_Kong2025_job
    submit_job_func = TMS_Kong2025.submit_TMS_Kong2025_job
    # Reduced-memory mode submits 1 job per unit of work (effectively uncapped);
    # otherwise jobs are capped at 50 to limit peak memory usage
    if reduce_mem_flag == "1":
        max_TMS_Kong2025_job = 99999
    else:
        max_TMS_Kong2025_job = 50
    try_count = 1
    try_max = 3  # inclusive of 1st try
    jobprefix = f"GTm_{pl_name}"

    # Extract `highest_split_bld` for `check_TMS_Kong2025_job` function
    gt_list_2_df = pd.read_csv(gt_list_2,
                               nrows=1,
                               dtype={
                                   'bld': str,
                                   'split_bld': str,
                                   'highest_split_bld': str
                               })
    highest_split_bld = gt_list_2_df['highest_split_bld'].iloc[0]

    check_kwargs = {
        'timestamp': timestamp,
        'try_count': try_count,
        'main_config_file': main_config_file,
        'pl_config_file': pl_config_file,
        'data_dir': data_dir,
        'dataset': dataset,
        'subID': subID,
        'run_type': run_type,
        'highest_new_bld': highest_split_bld,
        'target_type': target_type,
        'reduce_mem_flag': reduce_mem_flag
    }
    submit_kwargs = {
        'timestamp': timestamp,
        'jobprefix': jobprefix,
        'max_TMS_Kong2025_job': max_TMS_Kong2025_job,
        'main_config_file': main_config_file,
        'pl_config_file': pl_config_file,
        'dataset': dataset,
        'subID': subID,
        'run_type': run_type,
        'target_type': target_type,
        'reduce_mem_flag': reduce_mem_flag
    }

    # Check & re-submit job until hit try limit
    all_job_success = util.check_and_resubmit_job(check_job_func,
                                                  submit_job_func,
                                                  timestamp,
                                                  try_count,
                                                  try_max,
                                                  jobprefix,
                                                  username,
                                                  fail_list_flag=False,
                                                  check_kwargs=check_kwargs,
                                                  submit_kwargs=submit_kwargs)
