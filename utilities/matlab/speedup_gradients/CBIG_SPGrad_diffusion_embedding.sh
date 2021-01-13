#!/bin/bash
#############################################
# CBIG Speed Up Gradient Diffusion Emdedding 
#############################################
# AUTHOR #################################
# RU(BY) KONG 
# 2020/10/01  
##########################################
##########################################
# CBIG_SPGrad_diffusion_embedding.sh <output_dir> <num_component>
#
# In this script, we
# 1) Read in the gradient distance matrix: <output_dir>/?h_gradient_distance_matrix.npy
# 2) Perform diffusion embedding with <num_component> components on gradient distance matrix
# 3) Difffusion embeding results will be saved in output directory: 
#    <output_dir>/?h_emb_<num_component>_distance_matrix.mat
##########################################
# Example:
# CBIG_SPGrad_diffusion_embedding.sh /path/output_results 100
#
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

##########################
# Set up input arguments #
##########################
output_dir=$1;
num_component=$2;
lh_dist_mat="${output_dir}/lh_gradient_distance_matrix.npy";
rh_dist_mat="${output_dir}/rh_gradient_distance_matrix.npy";

if [ ! -f "$lh_dist_mat" ]; then
    echo "$lh_dist_mat does not exist!"
    exit 1
fi
if [ ! -f "$rh_dist_mat" ]; then
    echo "$rh_dist_mat does not exist!"
    exit 1
fi

dump_dir=$output_dir/dump
code_dir="$CBIG_CODE_DIR/utilities/matlab/speedup_gradients/utilities";

echo "Compute diffusion embedding matrix ..."
echo "NOTE: RuntimeWarning can be ignored"

if [ -f "$output_dir/rh_emb_${num_component}_distance_matrix.mat" ]; then
    echo "Embedding results already exist!"
else

    source activate CBIG_py3

    python $code_dir/apply_diffusion_on_individul_distance_matrix.py \
$lh_dist_mat $rh_dist_mat $output_dir $num_component

    if [ -f "$output_dir/rh_emb_${num_component}_distance_matrix.mat" ]; then
        echo "Finish computing diffusion embedding matrix!"
        touch $dump_dir/Done_diffusion_emdedding.dump
    else
        echo "Error: diffusion embedding matrix failed!"
    fi

    source deactivate CBIG_py3
fi