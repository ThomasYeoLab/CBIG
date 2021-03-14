#!/bin/bash
# this function runs replication of all results in our paper
#
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

root_dir=`dirname "$(readlink -f "$0")"`

bash root_dir/CBIG_pMFM_replication_part1_pMFM_main.sh
bash root_dir/CBIG_pMFM_replication_part2_pMFM_control_analysis.sh