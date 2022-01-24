#!/bin/bash
# this function generates all non-Python and non-R based figures for the TRBPC project
#
# Written by Angela Tam and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

outdir=$1
code_dir=$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/figure_utilities


# call matlab wrapper to generate all the PFMs
matlab -nodesktop -nosplash -nodisplay -r " addpath $code_dir; CBIG_TRBPC_replicate_matrix_plots( \
   '$outdir' ); exit; " 

# generate PFM similarity plots
source activate py_trbpc
in_dir=$code_dir/input
score_predicted=$in_dir/predicted_scores.txt
python $code_dir/clustering/multikernel_PFM_similarity_plots.py $in_dir $outdir $outdir $score_predicted
python $code_dir/clustering/similarity_fmri_clusters.py $outdir $outdir

# generate the cortical surface figures
declare -a cluster=("datadriven" "hypothesis" )
declare -a behavs=("cog" "ment_health" "pers" )
declare -a vals=("pos" "neg" )

for clust in ${cluster[@]}; do
for behav in ${behavs[@]}; do
for val in ${vals[@]}; do
datadir_fv="$outdir/$clust/conjunction/surfaces/${val}_${behav}"
echo $datadir_fv
$code_dir/surface_plots/CBIG_TRBPC_freeview_annot_fslr164k.sh $datadir_fv
done
done
done

# use R to generate the chord diagrams
source activate r_env

declare -a cluster=("datadriven" "hypothesis" )

for clust in ${cluster[@]}; do
datadir_r="$outdir/$clust/conjunction/chord"
$code_dir/chord/CBIG_TRBPC_chord.sh $datadir_r
#echo $datadir_r
done
