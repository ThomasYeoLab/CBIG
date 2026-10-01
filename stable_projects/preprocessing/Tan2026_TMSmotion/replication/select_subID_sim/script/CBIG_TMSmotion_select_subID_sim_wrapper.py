'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
"""
Examples:

python $TAN2026_TMSMOTION_DIR/replication/select_subID_sim/script/CBIG_TMSmotion_select_subID_sim_wrapper.py \
    $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json

python $TAN2026_TMSMOTION_DIR/replication/select_subID_sim/script/CBIG_TMSmotion_select_subID_sim_wrapper.py \
    $TAN2026_TMSMOTION_DIR/replication/pl2/config/pipeline.json

"""
# ===========================================================================
# IMPORTS
# ===========================================================================
import os
import sys
from datetime import datetime

TAN2026_TMSMOTION_DIR = os.environ.get('TAN2026_TMSMOTION_DIR')
sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'util'))
import CBIG_TMSmotion_util as util

sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'select_subID', 'sim'))
import CBIG_TMSmotion_select_subID_sim as sel_subID_sim

#############################################################################
# ===========================================================================
# GENERAL FOR ALL JOBS
# ===========================================================================
#############################################################################
# Read main config
main_config_file = os.path.join(TAN2026_TMSMOTION_DIR, 'replication', 'config',
                                'main_config.json')
main_config = util.read_json_file(main_config_file)
username = str(main_config['jobscheduler']['username'])

# Read pipeline config
pl_config_file = sys.argv[1]
pl_config = util.read_json_file(pl_config_file)
pl_name = str(pl_config['pipeline']['pl_name'])
rep_sel_subID_ngt_dir = pl_config['path']['rep_sel_subID_ngt_dir']
rep_sel_subID_sim_dir = pl_config['path']['rep_sel_subID_sim_dir']
subID_list = os.path.join(rep_sel_subID_ngt_dir, "subID_list_1.txt")

# Init convenience vars
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

#############################################################################
# ===========================================================================
# SIMULATE RUNS TO SELECT SUB_ID
# ===========================================================================
#############################################################################
# Init preproc job-specific vars
check_job_func = sel_subID_sim.check_all_select_subID_sim_job
submit_job_func = sel_subID_sim.submit_all_select_subID_sim_job
std_dir = os.path.join(rep_sel_subID_sim_dir, timestamp, 'job')

try_count = 1
try_max = 2  # inclusive of 1st try
jobprefix = f"SSS_{pl_name}"
check_kwargs = {
    'list_file': subID_list,
    'timestamp': timestamp,
    'try_count': try_count,
    'rep_sel_subID_sim_dir': rep_sel_subID_sim_dir
}
submit_kwargs = {
    'timestamp': timestamp,
    'jobprefix': jobprefix,
    'main_config_file': main_config_file,
    'pl_config_file': pl_config_file,
    'rep_sel_subID_sim_dir': rep_sel_subID_sim_dir,
    'list_file': subID_list,
    'std_dir': std_dir,
    'username': username
}
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
# GENERATE SELECTED SUB_ID SUMMARY
# ===========================================================================
#############################################################################
# Init preproc job-specific vars
check_job_func = sel_subID_sim.check_generate_select_subID_summary_job
submit_job_func = sel_subID_sim.submit_generate_select_subID_summary_job
std_dir = os.path.join(rep_sel_subID_sim_dir, timestamp, 'job')

try_count = 1
try_max = 2  # inclusive of 1st try
jobprefix = f"sum_{pl_name}"
check_kwargs = {
    'timestamp': timestamp,
    'try_count': try_count,
    'rep_sel_subID_sim_dir': rep_sel_subID_sim_dir
}
submit_kwargs = {
    'timestamp': timestamp,
    'jobprefix': jobprefix,
    'pl_config_file': pl_config_file,
    'rep_sel_subID_sim_dir': rep_sel_subID_sim_dir,
    'std_dir': std_dir,
    'username': username
}
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
