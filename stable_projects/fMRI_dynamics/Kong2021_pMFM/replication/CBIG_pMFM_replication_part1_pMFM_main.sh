#!/bin/bash
# this function runs replication of all results in our paper
#
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

root_dir=`dirname "$(readlink -f "$0")"`

cd root_dir/part1_pMFM_main/scripts
bash CBIG_pMFM_step1_training_main_wrapper.sh
bash CBIG_pMFM_step2_validation_main_wrapper.sh
bash CBIG_pMFM_step3_test_main_wrapper.sh
bash CBIG_pMFM_step4_generate_simulated_fc_fcd_main_wrapper.sh
bash CBIG_pMFM_step5_generate_STDFCD_correlation_main_wrapper.sh
bash CBIG_pMFM_step6_SWSTD_state_main_wrapper.sh
bash CBIG_pMFM_step7_perturbation_analysis_main_wrapper.sh
bash CBIG_pMFM_step8_gene_expression_analysis_desikan_wrapper.sh
bash CBIG_pMFM_step9_main_cleanup.sh
