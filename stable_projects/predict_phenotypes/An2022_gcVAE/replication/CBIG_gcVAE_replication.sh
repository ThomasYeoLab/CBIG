#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

## general setup
node_name=$(hostname)
if [ $node_name != 'headnode' ]; then
    echo "All replication jobs shoulde be submitted via headnode!"
    exit 1
fi
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
log_dir=$ROOTDIR'/job_logs'
mkdir -p $log_dir
job_prefix="GCVAE"
cd $ROOTDIR

## begin to replicate
echo ">>> Begin to reeplicate An2022_gcVAE..."
bash $ROOTDIR"/replication/scripts/CBIG_gcVAE_wait4jobs2finish.sh" $job_prefix

## step 0
echo ">>> Step0: Prepare for replication..."
bash $ROOTDIR"/replication/scripts/CBIG_gcVAE_replica_step0_prepare.sh"

## step 1
echo ">>> Step1: Matching datasets..."
s1_params=("ADNI-AIBL" "ADNI-MACC")
for s1_param_id in {0..1}; do
    s1_param=${s1_params[s1_param_id]}
    s1_name=${job_prefix}"_step1_"$s1_param
    s1_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step1_match_data.sh ${s1_param}"
    s1_err_log=${log_dir}"/"$s1_name".err"
    s1_out_log=${log_dir}"/"$s1_name".out"
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s1_cmd" \
        -name $s1_name -walltime 3:00:00 -mem 16G -ncpus 1 -joberr ${s1_err_log} -jobout ${s1_out_log}
    sleep 1s
done

## step 2
bash $ROOTDIR"/replication/scripts/CBIG_gcVAE_wait4jobs2finish.sh" $job_prefix
echo ">>> Step2: Split data into 10 folds..."
s2_name=${job_prefix}"_step2_split_data"
s2_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step2_split_data.sh"
s2_err_log=${log_dir}"/"$s2_name".err"
s2_out_log=${log_dir}"/"$s2_name".out"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s2_cmd" \
    -name $s2_name -walltime 2:00:00 -mem 8G -ncpus 1 -joberr ${s2_err_log} -jobout ${s2_out_log}
## step 2-1
bash $ROOTDIR"/replication/scripts/CBIG_gcVAE_wait4jobs2finish.sh" $job_prefix
s21_params=("ADNI-AIBL" "ADNI-MACC")
for s21_param_id in {0..1}; do
    s21_param=${s21_params[s21_param_id]}
    s21_name=${job_prefix}"_step2-1_"$s21_param
    s21_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step2-1_sample_data.sh ${s21_param}"
    s21_err_log=${log_dir}"/"$s21_name".err"
    s21_out_log=${log_dir}"/"$s21_name".out"
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s21_cmd" \
        -name $s21_name -walltime 2:00:00 -mem 16G -ncpus 1 -joberr ${s21_err_log} -jobout ${s21_out_log}
    sleep 1s
done

## setp 3
bash $ROOTDIR"/replication/scripts/CBIG_gcVAE_wait4jobs2finish.sh" $job_prefix
echo ">>> Step3: Generate input for unmatch2match & match2unmatch..."
s3_params_1=("ADNI-AIBL" "ADNI-MACC")
s3_params_2=("unmatch2match" "match2unmatch")
for s3_param_1_id in {0..1}; do
    for s3_param_2_id in {0..1}; do
        s3_param_1=${s3_params_1[s3_param_1_id]}
        s3_param_2=${s3_params_2[s3_param_2_id]}
        s3_name=${job_prefix}"_step3_"$s3_param_1"_"$s3_param_2
        s3_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step3_gen_input.sh ${s3_param_1} ${s3_param_2}"
        s3_err_log=${log_dir}"/"$s3_name".err"
        s3_out_log=${log_dir}"/"$s3_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s3_cmd" \
            -name $s3_name -walltime 2:00:00 -mem 8G -ncpus 1 -joberr ${s3_err_log} -jobout ${s3_out_log}
        sleep 1s
    done
done

## step 3-1
echo ">>> Step3-1: Generate input for sample_size..."
s31_params_1=("ADNI-AIBL" "ADNI-MACC")
s31_params_2=(10 20 30 40 50 60 70 80 90)
s31_params_3=(10 11 12 13 14 15 16 17 18 19)
for s31_param_1_id in {0..1}; do
    for s31_param_2_id in {0..8}; do
        for s31_param_3_id in {0..9}; do
            s31_param_1=${s31_params_1[s31_param_1_id]}
            s31_param_2=${s31_params_2[s31_param_2_id]}
            s31_param_3=${s31_params_3[s31_param_3_id]}
            s31_name=${job_prefix}"_step3-1_"$s31_param_1"_"$s31_param_2"_"$s31_param_3
            s31_cmd_1="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step3-1_gen_sampled_input.sh "
            s31_cmd=${s31_cmd_1}"${s31_param_1} ${s31_param_2} ${s31_param_3}"
            s31_err_log=${log_dir}"/"$s31_name".err"
            s31_out_log=${log_dir}"/"$s31_name".out"
            $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s31_cmd" \
                -name $s31_name -walltime 2:00:00 -mem 8G -ncpus 1 -joberr ${s31_err_log} -jobout ${s31_out_log}
            sleep 1s
        done
    done
done

## step 4
bash $ROOTDIR"/replication/scripts/CBIG_gcVAE_wait4jobs2finish.sh" $job_prefix
echo ">>> Step4: Train goalDNN model..."
## firstly train unmatch2match & match2unmatch experiments
s4_params_1=("ADNI-AIBL" "ADNI-MACC")
s4_params_2=("unmatch2match" "match2unmatch")
for s4_param_1_id in {0..1}; do
    for s4_param_2_id in {0..1}; do
        s4_param_1=${s4_params_1[s4_param_1_id]}
        s4_param_2=${s4_params_2[s4_param_2_id]}
        s4_name=${job_prefix}"_step4_"$s4_param_1"_"$s4_param_2
        s4_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step4_train_goalDNN.sh ${s4_param_1} ${s4_param_2}"
        s4_err_log=${log_dir}"/"$s4_name".err"
        s4_out_log=${log_dir}"/"$s4_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s4_cmd" \
            -name $s4_name -walltime 2:00:00 -mem 32G -ncpus 4 -ngpus 1 -joberr ${s4_err_log} -jobout ${s4_out_log}
        sleep 1s
    done
done
## then train sample_size experiments in parallel
s4_params_3=(10 20 30 40 50 60 70 80 90)
for s4_param_1_id in {0..1}; do
    s4_param_1=${s4_params_1[s4_param_1_id]}
    for s4_param_3_id in {0..8}; do
        # train models on GPU in parallel
        s4_param_3=${s4_params_3[s4_param_3_id]}
        s4_name=${job_prefix}"_step4-1_"$s4_param_1"_"$s4_param_3"perc"
        s4_cmd1="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step4-1_train_sampled_goalDNN.sh "
        s4_cmd=${s4_cmd1}"${s4_param_1} ${s4_param_3}"
        s4_err_log=${log_dir}"/"$s4_name".err"
        s4_out_log=${log_dir}"/"$s4_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s4_cmd" \
            -name $s4_name -walltime 12:00:00 -mem 32G -ncpus 4 -ngpus 1 -joberr ${s4_err_log} -jobout ${s4_out_log}
        sleep 1s
    done
done

# step 5
echo ">>> Step5: Fit ComBat model..."
## firstly train unmatch2match & match2unmatch experiments
s5_params_1=("ADNI-AIBL" "ADNI-MACC")
s5_params_2=("unmatch2match" "match2unmatch")
for s5_param_1_id in {0..1}; do
    for s5_param_2_id in {0..1}; do
        s5_param_1=${s5_params_1[s5_param_1_id]}
        s5_param_2=${s5_params_2[s5_param_2_id]}
        s5_name=${job_prefix}"_step5_"$s5_param_1"_"$s5_param_2
        s5_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step5_fit_ComBat.sh ${s5_param_1} ${s5_param_2}"
        s5_err_log=${log_dir}"/"$s5_name".err"
        s5_out_log=${log_dir}"/"$s5_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s5_cmd" \
            -name $s5_name -walltime 2:00:00 -mem 16G -ncpus 1 -joberr ${s5_err_log} -jobout ${s5_out_log}
        sleep 1s
    done
done
## then train sample_size experiments in parallel
s5_params_3=(10 20 30 40 50 60 70 80 90)
s5_params_4=(10 11 12 13 14 15 16 17 18 19)
for s5_param_1_id in {0..1}; do
    s5_param_1=${s5_params_1[s5_param_1_id]}
    for s5_param_3_id in {0..8}; do
        s5_param_3=${s5_params_3[s5_param_3_id]}
        for s5_param_4_id in {0..9}; do
            s5_param_4=${s5_params_4[s5_param_4_id]}
            s5_name=${job_prefix}"_step5-1_"$s5_param_1"_"$s5_param_3"perc_seed"$s5_param_4
            s5_cmd_1="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step5-1_fit_sampled_ComBat.sh "
            s5_cmd=${s5_cmd_1}"${s5_param_1} ${s5_param_3} ${s5_param_4}"
            s5_err_log=${log_dir}"/"$s5_name".err"
            s5_out_log=${log_dir}"/"$s5_name".out"
            $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s5_cmd" \
                -name $s5_name -walltime 2:00:00 -mem 16G -ncpus 1 -joberr ${s5_err_log} -jobout ${s5_out_log}
            sleep 1s
        done
    done
done

## step 6
echo ">>> Step6: Train cVAE model..."
## firstly train unmatch2match & match2unmatch experiments
s6_params_1=("ADNI-AIBL" "ADNI-MACC")
s6_params_2=("unmatch2match" "match2unmatch")
for s6_param_1_id in {0..1}; do
    for s6_param_2_id in {0..1}; do
        s6_param_1=${s6_params_1[s6_param_1_id]}
        s6_param_2=${s6_params_2[s6_param_2_id]}
        s6_name=${job_prefix}"_step6_"$s6_param_1"_"$s6_param_2
        s6_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step6_train_cVAE.sh ${s6_param_1} ${s6_param_2}"
        s6_err_log=${log_dir}"/"$s6_name".err"
        s6_out_log=${log_dir}"/"$s6_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s6_cmd" \
            -name $s6_name -walltime 2:00:00 -mem 32G -ncpus 4 -ngpus 1 -joberr ${s6_err_log} -jobout ${s6_out_log}
        sleep 1s
    done
done
## then train sample_size experiments in parallel
s6_params_3=(10 20 30 40 50 60 70 80 90)
for s6_param_1_id in {0..1}; do
    s6_param_1=${s6_params_1[s6_param_1_id]}
    for s6_param_3_id in {0..8}; do
        # train models on GPU in parallel
        s6_param_3=${s6_params_3[s6_param_3_id]}
        s6_name=${job_prefix}"_step6-1_"$s6_param_1"_"$s6_param_3"perc"
        s6_cmd_1="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step6-1_train_sampled_cVAE.sh "
        s6_cmd=${s6_cmd_1}"${s6_param_1} ${s6_param_3}"
        s6_err_log=${log_dir}"/"$s6_name".err"
        s6_out_log=${log_dir}"/"$s6_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s6_cmd" \
            -name $s6_name -walltime 12:00:00 -mem 32G -ncpus 4 -ngpus 1 -joberr ${s6_err_log} -jobout ${s6_out_log}
        sleep 1s
    done
done

## step 7
bash $ROOTDIR"/replication/scripts/CBIG_gcVAE_wait4jobs2finish.sh" $job_prefix
echo ">>> Step7: Train gcVAE model..."
## firstly train unmatch2match & match2unmatch experiments
s7_params_1=("ADNI-AIBL" "ADNI-MACC")
s7_params_2=("unmatch2match" "match2unmatch")
for s7_param_1_id in {0..1}; do
    for s7_param_2_id in {0..1}; do
        s7_param_1=${s7_params_1[s7_param_1_id]}
        s7_param_2=${s7_params_2[s7_param_2_id]}
        s7_name=${job_prefix}"_step7_"$s7_param_1"_"$s7_param_2
        s7_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step7_train_gcVAE.sh ${s7_param_1} ${s7_param_2}"
        s7_err_log=${log_dir}"/"$s7_name".err"
        s7_out_log=${log_dir}"/"$s7_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s7_cmd" \
            -name $s7_name -walltime 2:00:00 -mem 32G -ncpus 4 -ngpus 1 -joberr ${s7_err_log} -jobout ${s7_out_log}
        sleep 1s
    done
done
## then train sample_size experiments in parallel
s7_params_3=(10 20 30 40 50 60 70 80 90)
for s7_param_1_id in {0..1}; do
    s7_param_1=${s7_params_1[s7_param_1_id]}
    for s7_param_3_id in {0..8}; do
        # train models on GPU in parallel
        s7_param_3=${s7_params_3[s7_param_3_id]}
        s7_name=${job_prefix}"_step7-1_"$s7_param_1"_"$s7_param_3"perc"
        s7_cmd_1="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step7-1_train_sampled_gcVAE.sh "
        s7_cmd=${s7_cmd_1}"${s7_param_1} ${s7_param_3}"
        s7_err_log=${log_dir}"/"$s7_name".err"
        s7_out_log=${log_dir}"/"$s7_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s7_cmd" \
            -name $s7_name -walltime 12:00:00 -mem 32G -ncpus 4 -ngpus 1 -joberr ${s7_err_log} -jobout ${s7_out_log}
        sleep 1s
    done
done

## step 8
bash $ROOTDIR"/replication/scripts/CBIG_gcVAE_wait4jobs2finish.sh" $job_prefix
echo ">>> Step8: Harmonize data using trained VAE models..."
## firstly train unmatch2match & match2unmatch experiments
s8_params_1=("ADNI-AIBL" "ADNI-MACC")
s8_params_2=("unmatch2match" "match2unmatch")
for s8_param_1_id in {0..1}; do
    for s8_param_2_id in {0..1}; do
        s8_param_1=${s8_params_1[s8_param_1_id]}
        s8_param_2=${s8_params_2[s8_param_2_id]}
        s8_name=${job_prefix}"_step8_"$s8_param_1"_"$s8_param_2
        s8_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step8_perform_VAE_harm.sh ${s8_param_1} ${s8_param_2}"
        s8_err_log=${log_dir}"/"$s8_name".err"
        s8_out_log=${log_dir}"/"$s8_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s8_cmd" \
            -name $s8_name -walltime 4:00:00 -mem 32G -ncpus 1 -joberr ${s8_err_log} -jobout ${s8_out_log}
        sleep 1s
    done
done
## then train sample_size experiments in parallel
s8_params_3=(10 20 30 40 50 60 70 80 90)
s8_params_4=(10 11 12 13 14 15 16 17 18 19)
for s8_param_1_id in {0..1}; do
    s8_param_1=${s8_params_1[s8_param_1_id]}
    for s8_param_3_id in {0..8}; do
        s8_param_3=${s8_params_3[s8_param_3_id]}
        for s8_param_4_id in {0..9}; do
            s8_param_4=${s8_params_4[s8_param_4_id]}
            s8_name=${job_prefix}"_step8-1_"$s8_param_1"_"$s8_param_3"perc_seed"$s8_param_4
            s8_cmd_1="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step8-1_perform_sampled_VAE_harm.sh "
            s8_cmd=${s8_cmd_1}"${s8_param_1} ${s8_param_3} ${s8_param_4}"
            s8_err_log=${log_dir}"/"$s8_name".err"
            s8_out_log=${log_dir}"/"$s8_name".out"
            $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s8_cmd" \
                -name $s8_name -walltime 4:00:00 -mem 32G -ncpus 1 -joberr ${s8_err_log} -jobout ${s8_out_log}
            sleep 1s
        done
    done
done

## step 9
bash $ROOTDIR"/replication/scripts/CBIG_gcVAE_wait4jobs2finish.sh" $job_prefix
echo ">>> Step9:Train XGBoost model for dataset prediction..."
## firstly train unmatch2match & match2unmatch experiments
s9_params_1=("ADNI-AIBL" "ADNI-MACC")
s9_params_2=("unmatch2match" "match2unmatch")
for s9_param_1_id in {0..1}; do
    for s9_param_2_id in {0..1}; do
        s9_param_1=${s9_params_1[s9_param_1_id]}
        s9_param_2=${s9_params_2[s9_param_2_id]}
        s9_name=${job_prefix}"_step9_"$s9_param_1"_"$s9_param_2
        s9_cmd="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step9_train_XGBoost.sh ${s9_param_1} ${s9_param_2}"
        s9_err_log=${log_dir}"/"$s9_name".err"
        s9_out_log=${log_dir}"/"$s9_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s9_cmd" \
            -name $s9_name -walltime 6:00:00 -mem 16G -ncpus 1 -joberr ${s9_err_log} -jobout ${s9_out_log}
        sleep 1s
    done
done
## then train sample_size experiments in parallel
s9_params_3=(10 20 30 40 50 60 70 80 90)
s9_params_4=(10 11 12 13 14 15 16 17 18 19)
for s9_param_1_id in {0..1}; do
    s9_param_1=${s9_params_1[s9_param_1_id]}
    for s9_param_3_id in {0..8}; do
        s9_param_3=${s9_params_3[s9_param_3_id]}
        for s9_param_4_id in {0..9}; do
            s9_param_4=${s9_params_4[s9_param_4_id]}
            s9_name=${job_prefix}"_step9-1_"$s9_param_1"_"$s9_param_3"perc_seed"$s9_param_4
            s9_cmd_1="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step9-1_train_sampled_XGBoost.sh "
            s9_cmd=${s9_cmd_1}"${s9_param_1} ${s9_param_3} ${s9_param_4}"
            s9_err_log=${log_dir}"/"$s9_name".err"
            s9_out_log=${log_dir}"/"$s9_name".out"
            $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s9_cmd" \
                -name $s9_name -walltime 6:00:00 -mem 16G -ncpus 1 -joberr ${s9_err_log} -jobout ${s9_out_log}
            sleep 1s
        done
    done
done

## step 10
echo ">>> Step10:Evaluate downstream application performance after harmonization using goalDNN..."
## firstly train unmatch2match & match2unmatch experiments
s10_params_1=("ADNI-AIBL" "ADNI-MACC")
s10_params_2=("unmatch2match" "match2unmatch")
for s10_param_1_id in {0..1}; do
    for s10_param_2_id in {0..1}; do
        s10_param_1=${s10_params_1[s10_param_1_id]}
        s10_param_2=${s10_params_2[s10_param_2_id]}
        s10_name=${job_prefix}"_step10_"$s10_param_1"_"$s10_param_2
        s10_cmd_1="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step10_goalDNN_evaluation.sh "
        s10_cmd=${s10_cmd_1}"${s10_param_1} ${s10_param_2}"
        s10_err_log=${log_dir}"/"$s10_name".err"
        s10_out_log=${log_dir}"/"$s10_name".out"
        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s10_cmd" \
            -name $s10_name -walltime 2:00:00 -mem 16G -ncpus 1 -joberr ${s10_err_log} -jobout ${s10_out_log}
        sleep 1s
    done
done
## then train sample_size experiments in parallel
s10_params_3=(10 20 30 40 50 60 70 80 90)
s10_params_4=(10 11 12 13 14 15 16 17 18 19)
for s10_param_1_id in {0..1}; do
    s10_param_1=${s10_params_1[s10_param_1_id]}
    for s10_param_3_id in {0..8}; do
        s10_param_3=${s10_params_3[s10_param_3_id]}
        for s10_param_4_id in {0..9}; do
            s10_param_4=${s10_params_4[s10_param_4_id]}
            s10_name=${job_prefix}"_step10-1_"$s10_param_1"_"$s10_param_3"perc_seed"$s10_param_4
            s10_cmd_1="$ROOTDIR/replication/scripts/CBIG_gcVAE_replica_step10-1_goalDNN_sampled_evaluation.sh "
            s10_cmd=${s10_cmd_1}"${s10_param_1} ${s10_param_3} ${s10_param_4}"
            s10_err_log=${log_dir}"/"$s10_name".err"
            s10_out_log=${log_dir}"/"$s10_name".out"
            $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$s10_cmd" \
                -name $s10_name -walltime 2:00:00 -mem 16G -ncpus 1 -joberr ${s10_err_log} -jobout ${s10_out_log}
            sleep 1s
        done
    done
done
