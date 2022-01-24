# this function will generate the chord diagrams for the TRBPC project

# Written by Nanbo Sun, Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# path to input and outputs
datadir=$1

code_dir=$CBIG_CODE_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/figure_utilities/chord

declare -a clus=("cog" "ment_health" "pers" )
declare -a vals=("abs" "pos" "neg" "bothdir" )

for clu in ${clus[@]}; do
for val in ${vals[@]}; do
input="$datadir/avg_${val}_${clu}.csv"
link_colorscale=bwr
min_thre=0.00000001
max_thre=1.0
output="$datadir/chord_${val}_${clu}_${min_thre}_${max_thre}"

R CMD BATCH --no-save --no-restore "--args $input $link_colorscale $min_thre $max_thre $output" \
$code_dir/CBIG_TRBPC_chord_diagram_18x18.r chord.Rout
done
done
