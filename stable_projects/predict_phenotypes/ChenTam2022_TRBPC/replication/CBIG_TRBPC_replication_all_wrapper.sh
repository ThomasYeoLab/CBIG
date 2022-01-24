#!/bin/bash
# this function runs replication of all results in our paper
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

outdir=$1
root_dir=`dirname "$(readlink -f "$0")"`

# replicate regression models of main analysis
$root_dir/CBIG_TRBPC_regression_main_wrapper.sh $outdir/KRR

# replicate regression models of control analysis
$root_dir/CBIG_TRBPC_regression_control_wrapper.sh $outdir

# replicate predictive feature matrix of all models
$root_dir/CBIG_TRBPC_PFM_wrapper.sh $outdir $outdir/PFM

# replicate permutation test of KRR
$root_dir/CBIG_TRBPC_KRR_perm_test_wrapper.sh $outdir/KRR $outdir/perm_test/KRR

# replicate permutation test of PFMS
$root_dir/CBIG_TRBPC_PFM_perm_test_wrapper.sh $outdir/KRR/allFC $outdir/perm_test/PFM