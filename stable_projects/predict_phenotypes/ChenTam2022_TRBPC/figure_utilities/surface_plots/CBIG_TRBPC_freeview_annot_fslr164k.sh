# this function plots values on the cortical surface with freeview

# Written by Ruby Kong, Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

data_dir=$1
id=$2
annot_dir=$data_dir/annot
fig_dir=$data_dir/figures
mkdir -p $fig_dir
cd $annot_dir

for hemi in {lh,rh}; do

for lh in ${hemi}*$id.annot; do
	if [ $hemi == "lh" ]; then
		freeview -f \
		$CBIG_CODE_DIR/data/templates/surface/fs_LR_164k/surf/${hemi}.very_inflated:annot=$annot_dir/$lh:edgethickness=0 \
		-zoom 1.78 -ss $fig_dir/$lh.lateral.png
	else
		freeview -f \
		$CBIG_CODE_DIR/data/templates/surface/fs_LR_164k/surf/${hemi}.very_inflated:annot=$annot_dir/$lh:edgethickness=0 \
		-cam Azimuth 180 -zoom 1.78 -ss $fig_dir/$lh.lateral.png
	fi
	if [ $hemi == "lh" ]; then
		freeview -f \
		$CBIG_CODE_DIR/data/templates/surface/fs_LR_164k/surf/${hemi}.very_inflated:annot=$annot_dir/$lh:edgethickness=0 \
		-cam Azimuth 180 -zoom 1.68 -ss $fig_dir/$lh.medial.png
	else
		freeview -f \
		$CBIG_CODE_DIR/data/templates/surface/fs_LR_164k/surf/${hemi}.very_inflated:annot=$annot_dir/$lh:edgethickness=0 \
		-zoom 1.68 -ss $fig_dir/$lh.medial.png
	fi
	
done
done
