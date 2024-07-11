#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

## general setup
node_name=$(hostname)
if [ $node_name != 'headnode' ]; then
	echo "All replication jobs shoulde be submitted via headnode!"
	exit 1
fi
base_dir=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/Chen2024_MMM"
script_dir=$base_dir"/scripts"
rep_dir=$base_dir"/replication"
log_dir=$rep_dir"/log"
mkdir -p $log_dir
job_prefix="MMM"
cd $rep_dir

echo ">>> Begin to replicate Chen2024_MMM..."
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix

## step 1: run classical KRR on target dataset (HCP-YA and HCP-Aging)
echo ">>> Step1: classical KRR ..."
s1_params=("HCP" "HCPA")
for tar_dataset in ${s1_params[@]}; do
	sh $script_dir"/CBIG_MMM_KRR_classical_wrapper.sh" ${tar_dataset}
	sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix
	# summarize classical KRR results
	sh $script_dir"/CBIG_MMM_KRR_classical_summary.sh" ${tar_dataset}
done

## step 2: get data split files for target datasets (HCP-YA, HCP-Aging)
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix
echo ">>> Step 2: get data split files ..."
s2_params=("HCP" "HCPA")
for tar_dataset in ${s2_params[@]}; do
	sh $script_dir"/CBIG_MMM_get_split.sh" ${tar_dataset}
done

## step 3: train DNN base model on extra-large dataset (UK Biobank) and predict on target datasets (HCP-YA, HCP-Aging)
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix
echo ">>> Step 3: DNN base model ..."
# Train DNN model on UK Biobank. The hyperparameters are automaticaly tuned by optuna
sh $script_dir"/CBIG_MMM_DNN_xlarge_train_submit_job.sh" "UKBB"
# Use DNN models (trained on UK Biobank) to predict phenotypes for target datasets
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix"_DNN_train"
sh $script_dir"/CBIG_MMM_DNN_xlarge_predict_submit_job.sh" "UKBB"

## step 4: fine-tune DNN model (pre-trained on UK Biobank) to predict target phenotypes on HCP-YA & HCP-Aging
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix"_DNN_train"
echo ">>> Step 4: transfer learning to target datasets"
s4_params=("HCP" "HCPA")
for tar_dataset in ${s4_params[@]}; do
	sh $script_dir"/CBIG_MMM_transfer_learning_submit_job.sh" "UKBB" ${tar_dataset}
	sleep 1s
done

## step 5: train LRR base model on all source dataset and predict on target dataset (HCP-YA, HCP-Aging)
echo ">>> Step 5: LRR base model (layer 1)..."
sh $script_dir"/CBIG_MMM_RR_xlarge_train_submit_jobs.sh" "UKBB"
sh $script_dir"/CBIG_MMM_RR_large_submit_jobs.sh" "ABCD" "1layer"
sh $script_dir"/CBIG_MMM_RR_medium_submit_jobs.sh" "GSP" "1layer"
sh $script_dir"/CBIG_MMM_RR_medium_submit_jobs.sh" "HBN" "1layer"
sh $script_dir"/CBIG_MMM_RR_medium_submit_jobs.sh" "eNKI" "1layer"
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix"_RR_train_XL"
sh $script_dir"/CBIG_MMM_RR_xlarge_predict_submit_job.sh" "UKBB"

## step 6: train KRR stacking model on large & medium source datasets and predict on target dataset (HCP-YA, HCP-Aging)
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix"_DNN_predict"
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix"_RR"
echo ">>> Step 6: KRR stacking model on source datasets (layer 2) ..."
# Apply meta-matching from extra-large source dataset (UK Biobank),
# then train KRR on large source dataset (ABCD) and apply on HCP-YA & HCP-Aging target datasets
sh $script_dir"/CBIG_MMM_RR_large_submit_jobs.sh" "ABCD" "2layer"
# Apply meta-matching from extra-large source dataset (UK Biobank) + large source dataset (ABCD),
# then train KRR on medium source datasets (GSP, HBN, eNKI) and apply on HCP-YA & HCP-Aging target datasets
sh $script_dir"/CBIG_MMM_RR_medium_submit_jobs.sh" "GSP" "2layer"
sh $script_dir"/CBIG_MMM_RR_medium_submit_jobs.sh" "HBN" "2layer"
sh $script_dir"/CBIG_MMM_RR_medium_submit_jobs.sh" "eNKI" "2layer"

## step 7: run meta-matching methods
# (meta-matching with stacking, meta-matching with dataset stacking, multilayer meta-matching)
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix"_DNN"
sh $script_dir"/CBIG_MMM_wait4jobs2finish.sh" $job_prefix"_RR"
echo ">>> Step 7: Run meta-matching methods ..."
s7_datasets=("HCP" "HCPA")
## Meta-matching with stacking, meta-matching with dataset stacking, multilayer meta-matching
s7_methods=("MM_stacking" "dataset_stacking" "multilayer_stacking")
for dataset in ${s7_datasets[@]}; do
	for method in ${s7_methods[@]}; do
		sh $script_dir"/CBIG_MMM_stacking_submit_job.sh" $method $dataset
	done
done
