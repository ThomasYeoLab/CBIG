#!/bin/bash

# Replication wrapper script
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

## general setup
node_name=$(hostname)
if [ $node_name != 'headnode' ]; then
    echo "Error: All replication jobs should be submitted via headnode!"
    exit 1
fi
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
log_dir=$ROOTDIR'/job_logs'
mkdir -p $log_dir
job_prefix="DRB"
cd $ROOTDIR

## replication starts
echo ">>> Begin to reeplicate An2024_DeepResBat..."
bash $ROOTDIR"/replication/scripts/CBIG_DeepResBat_wait4job2finish.sh" $job_prefix

################################### step 0 ####################################
echo ">>> Step0: Prepare for replication..."
bash $ROOTDIR"/replication/scripts/CBIG_DeepResBat_replica_step0_prepare.sh"
bash $ROOTDIR"/replication/scripts/CBIG_DeepResBat_wait4job2finish.sh" $job_prefix
################################### step 0 ####################################

################################### step 1 ####################################
echo ">>> Step1: Generate input data..."
s1_cmd="bash $ROOTDIR/replication/scripts/CBIG_DeepResBat_replica_step1_gen_input.sh ADNI-AIBL"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s1_cmd" \
    -name $job_prefix"_Step1_ADNI-AIBL" \
    -walltime 02:00:00 \
    -mem 2G \
    -ncpus 1 \
    -jobout $log_dir"/DRB_Step1_ADNI-AIBL.out" \
    -joberr $log_dir"/DRB_Step1_ADNI-AIBL.err"
sleep 5m
s1_cmd="bash $ROOTDIR/replication/scripts/CBIG_DeepResBat_replica_step1_gen_input.sh ADNI-MACC"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s1_cmd" \
    -name $job_prefix"_Step1_ADNI-MACC" \
    -walltime 02:00:00 \
    -mem 2G \
    -ncpus 1 \
    -jobout $log_dir"/DRB_Step1_ADNI-MACC.out" \
    -joberr $log_dir"/DRB_Step1_ADNI-MACC.err"
bash $ROOTDIR"/replication/scripts/CBIG_DeepResBat_wait4job2finish.sh" $job_prefix
################################### step 1 ####################################

################################### step 2 ####################################
echo ">>> Step2: Baseline harmonization models: ComBat & CovBat...."
s2_cmd="bash $ROOTDIR/replication/scripts/CBIG_DeepResBat_replica_step2_fit_ComBats.sh ADNI-AIBL"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s2_cmd" \
    -name $job_prefix"_Step2_ADNI-AIBL" \
    -walltime 02:00:00 \
    -mem 2G \
    -ncpus 1 \
    -jobout $log_dir"/DRB_Step2_ADNI-AIBL.out" \
    -joberr $log_dir"/DRB_Step2_ADNI-AIBL.err"
sleep 5s
s2_cmd="bash $ROOTDIR/replication/scripts/CBIG_DeepResBat_replica_step2_fit_ComBats.sh ADNI-MACC"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s2_cmd" \
    -name $job_prefix"_Step2_ADNI-MACC" \
    -walltime 02:00:00 \
    -mem 2G \
    -ncpus 1 \
    -jobout $log_dir"/DRB_Step2_ADNI-MACC.out" \
    -joberr $log_dir"/DRB_Step2_ADNI-MACC.err"
################################### step 2 ####################################

################################### step 3 ####################################
echo ">>> Step3: Baseline harmonization models: cVAE...."
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step3.sh "ADNI-AIBL"
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step3.sh "ADNI-MACC"
################################### step 3 ####################################

################################### step 4 ####################################
echo ">>> Step4: Proposed harmonization models: coVAE...."
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step4.sh "ADNI-AIBL"
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step4.sh "ADNI-MACC"
################################### step 4 ####################################

################################### step 5 ####################################
echo ">>> Step5: Proposed harmonization models: DeepResBat...."
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step5_ADNI-AIBL.sh
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step5_ADNI-MACC.sh
################################### step 5 ####################################

bash $ROOTDIR"/replication/scripts/CBIG_DeepResBat_wait4job2finish.sh" $job_prefix

################################### step 6 ####################################
echo ">>> Step6: Evaluate performance by predicting dataset...."
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step6.sh "ADNI-AIBL"
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step6.sh "ADNI-MACC"
################################### step 6 ####################################

################################### step 7 ####################################
echo ">>> Step7: Evaluate performance by MANOVA...."
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step7.sh "ADNI-AIBL"
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step7.sh "ADNI-MACC"
################################### step 7 ####################################

################################### step 8 ####################################
echo ">>> Step8: Evaluate performance by GLM...."
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step8.sh "ADNI-AIBL"
bash $ROOTDIR/replication/submission/CBIG_DeepResBat_step8.sh "ADNI-MACC"
################################### step 8 ####################################
