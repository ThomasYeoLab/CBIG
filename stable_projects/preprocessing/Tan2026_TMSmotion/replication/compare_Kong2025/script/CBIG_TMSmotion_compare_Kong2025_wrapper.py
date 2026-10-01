'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
"""
Examples:

python $TAN2026_TMSMOTION_DIR/replication/compare_Kong2025/script/CBIG_TMSmotion_compare_Kong2025_wrapper.py \
    $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json

python $TAN2026_TMSMOTION_DIR/replication/compare_Kong2025/script/CBIG_TMSmotion_compare_Kong2025_wrapper.py \
    $TAN2026_TMSMOTION_DIR/replication/pl2/config/pipeline.json

"""
#############################################################################
# ===========================================================================
# IMPORTS
# ===========================================================================
#############################################################################

import os
import sys
from datetime import datetime

TAN2026_TMSMOTION_DIR = os.environ.get('TAN2026_TMSMOTION_DIR')
sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'util'))
import CBIG_TMSmotion_util as util

CBIG_CODE_DIR = util.get_cbig_code_dir()

sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'compare'))
import CBIG_TMSmotion_compare_Kong2025 as compare_Kong2025

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
rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']

# Define grid or slide mode
mode = 'grid'

# Define TMS type
TMS_type = 'Kong2025_dep'

# Define input files
list_file = os.path.join(rep_sel_subID_sim_dir, "result", "summary",
                         "selected_subID_list_2.txt")

# Init convenience vars
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

# Init low/high memory pipeline
reduce_mem_flag = "1"

# Init network type
network_type = "17_network"

#############################################################################
# ===========================================================================
# TMS TARGET EUCLIDEAN DISTANCE
# ===========================================================================
#############################################################################
check_job_func = compare_Kong2025.check_compare_Kong2025_TMS_wrapper_job
submit_job_func = compare_Kong2025.submit_compare_Kong2025_TMS_wrapper_job
try_count = 1
try_max = 1  # inclusive of 1st try
jobprefix = f"CTm_{pl_name}"

check_kwargs = {
    'timestamp': timestamp,
    'try_count': try_count,
    'main_config_file': main_config_file,
    'pl_config_file': pl_config_file,
    'TMS_type': TMS_type
}
submit_kwargs = {
    'timestamp': timestamp,
    'jobprefix': jobprefix,
    'main_config_file': main_config_file,
    'pl_config_file': pl_config_file,
    'list_file': list_file,
    'mode': mode,
    'TMS_type': TMS_type
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
# MSHBM DICE COEFFICIENT
# ===========================================================================
#############################################################################
check_job_func = compare_Kong2025.check_compare_Kong2025_MSHBM_wrapper_job
submit_job_func = compare_Kong2025.submit_compare_Kong2025_MSHBM_wrapper_job
try_count = 1
try_max = 1  # inclusive of 1st try
jobprefix = f"CMs_{pl_name}"

check_kwargs = {
    'timestamp': timestamp,
    'try_count': try_count,
    'main_config_file': main_config_file,
    'pl_config_file': pl_config_file
}
submit_kwargs = {
    'timestamp': timestamp,
    'jobprefix': jobprefix,
    'main_config_file': main_config_file,
    'pl_config_file': pl_config_file,
    'list_file': list_file,
    'network_type': network_type,
    'mode': mode
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
# SUMMARIZE COMPARISON OF OUTPUT TYPES IN SINGLE DATAFRAME
# ===========================================================================
#############################################################################
check_job_func = compare_Kong2025.check_all_generate_compare_summary_subID_list_file_job
submit_job_func = compare_Kong2025.submit_all_generate_compare_summary_subID_list_file_job
try_count = 1
try_max = 1  # inclusive of 1st try
jobprefix = f"CSu_{pl_name}"

check_kwargs = {
    'timestamp': timestamp,
    'try_count': try_count,
    'pl_config_file': pl_config_file,
    'list_file': list_file,
    'TMS_type': TMS_type
}
submit_kwargs = {
    'timestamp': timestamp,
    'jobprefix': jobprefix,
    'pl_config_file': pl_config_file,
    'list_file': list_file,
    'mode': mode,
    'TMS_type': TMS_type
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
# CONCATENATE COMPARISON SUMMARY ACROSS SUBIDS
# ===========================================================================
#############################################################################
check_job_func = compare_Kong2025.check_concat_all_compare_summary_subID_list_file_job
submit_job_func = compare_Kong2025.submit_concat_all_compare_summary_subID_list_file_job
try_count = 1
try_max = 1  # inclusive of 1st try
jobprefix = f"CCc_{pl_name}"

check_kwargs = {
    'timestamp': timestamp,
    'try_count': try_count,
    'pl_config_file': pl_config_file,
    'TMS_type': TMS_type
}
submit_kwargs = {
    'timestamp': timestamp,
    'jobprefix': jobprefix,
    'pl_config_file': pl_config_file,
    'list_file': list_file,
    'TMS_type': TMS_type
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
