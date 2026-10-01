'''
Written by Trevor Tan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
"""
Examples:

python $TAN2026_TMSMOTION_DIR/replication/select_subID_window/script/CBIG_TMSmotion_select_subID_window_wrapper.py \
    $TAN2026_TMSMOTION_DIR/replication/pl1/config/pipeline.json

python $TAN2026_TMSMOTION_DIR/replication/select_subID_window/script/CBIG_TMSmotion_select_subID_window_wrapper.py \
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

sys.path.append(os.path.join(TAN2026_TMSMOTION_DIR, 'select_subID', 'window'))
import CBIG_TMSmotion_select_subID_window as sel_subID_window

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
rep_sel_subID_window_dir = pl_config['path']['rep_sel_subID_window_dir']

# Init convenience vars
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

#############################################################################
# ===========================================================================
# GENERATE WINDOW LIST
# ===========================================================================
#############################################################################
check_job_func = sel_subID_window.check_select_subID_window_job
submit_job_func = sel_subID_window.submit_select_subID_window_job
std_dir = os.path.join(rep_sel_subID_window_dir, timestamp, 'job')

try_count = 1
try_max = 1  # inclusive of 1st try
jobprefix = f"SWn_{pl_name}"

check_kwargs = {
    'timestamp': timestamp,
    'try_count': try_count,
    'rep_sel_subID_window_dir': rep_sel_subID_window_dir
}
submit_kwargs = {
    'timestamp': timestamp,
    'jobprefix': jobprefix,
    'main_config_file': main_config_file,
    'pl_config_file': pl_config_file,
    'rep_sel_subID_ngt_dir': rep_sel_subID_ngt_dir,
    'rep_sel_subID_window_dir': rep_sel_subID_window_dir,
    'std_dir': std_dir
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
