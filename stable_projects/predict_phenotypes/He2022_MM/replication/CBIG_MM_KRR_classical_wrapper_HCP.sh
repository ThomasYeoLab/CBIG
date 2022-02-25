#!/bin/sh
# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Initialize directory
input_dir="$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/He2022_MM"
base_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM"
code_dir="${base_dir}/KRR_CLASSICAL"
rep_dir="${base_dir}/replication"

# Initialize input and output
ks=(10 20 50 100 200)
split="diff_roi"
final_phe_list="${input_dir}/HCP_${split}_final_phe_list.txt"
phe_csv="${input_dir}/HCP_${split}_final.csv"
subj_list="${input_dir}/HCP_${split}_subj_list.txt"
fc_mat="${input_dir}/HCP_${split}_pfc.mat"
output="${rep_dir}/output_KRR_classical_HCP"

# Run KRR by submitting job
for k in "${ks[@]}"; do
	while read p; do
		n_job="$(qstat -u the | wc -l)"
		while [ $n_job -gt 50 ]; do
			echo $n_job large than 50, sleep for 1 mins
			sleep 1m
			n_job="$(qstat -u the | wc -l)"
		done
		sh CBIG_MM_KRR_classical_submit_job.sh $p $output $phe_csv $code_dir $k $split $subj_list $fc_mat
	done <$final_phe_list
done
