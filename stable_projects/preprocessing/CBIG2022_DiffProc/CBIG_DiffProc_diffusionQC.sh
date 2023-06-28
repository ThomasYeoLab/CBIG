#!/bin/sh
#####
# Example: 
#    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_diffusion_processing2022/CBIG_diffusionQC.sh \
#        --subj $subject --subjDIR $projdir/$subject/T1w/Diffusion \
#        --subjImg data.nii.gz --outdir $base_outdir/S1200/individuals/$subject \
#        --bval_file bvals --bvec_file bvecs \
#        --brainmask $projdir/$subject/T1w/Diffusion/nodif_brain_mask.nii.gz \
#        --round_bvals_1000 --use_single_shell 1000 --del_files
#
# This function runs QC for diffusion data and fits DTItensor to specified subject. Results are saved in a txt file.
# Intervolume motion, the sum of squared differences of each shell, b0 SNR, and residual fit to DTI model. 
#
# Written by Leon Ooi.
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

#####
# parse args
#####
scriptdir=$( dirname "${BASH_SOURCE[0]}" )
while [[ $# -gt 0 ]]; do
key="$1"

case $key in

    -s|--subj)
    subj="$2"
        shift; shift;;

    -d|--subjDIR)
    subjDIR="$2"
        shift; shift;;

    -i|--subjImg)
    subjImg="$2"
        shift; shift;;

    -a|--bval_file)
    bval_file="$2"
        shift; shift;;

    -e|--bvec_file)
    bvec_file="$2"
        shift; shift;;

    -o|--outdir)
    outdir="$2"
        shift; shift;;

    -b|--brainmask)
    brainmask="$2"
        shift; shift;;

    -r|--round_bvals_1000)
    round_bvals_1000="y"
        shift;;

    -u|--use_single_shell)
    use_single_shell="$2"
        shift; shift;;

    -d|--del_files)
    del_files="y"
        shift;;

    *)    # unknown option
    echo "Unknown option: $1"
    shift;;
esac
done

# check args
if [ -z "$subj" ] || [ -z "$subjDIR" ] || [ -z "$subjImg" ] \
|| [ -z "$bval_file" ] || [ -z "$bvec_file" ] || [ -z "$outdir" ]; then
    echo "Missing compulsory variable!"
    exit
fi

if [ -z "$brainmask" ]; then
    brainmask="NIL"
fi

if [ -z "$round_bvals_1000" ]; then
    round_bvals_1000="n"
fi

if [ -z "$use_single_shell" ]; then
    use_single_shell="NIL"
fi

if [ -z "$del_files" ]; then
    del_files="n"
fi

# print selected options
echo "[subj] = $subj"
echo "[subjDIR] = $subjDIR"
echo "[subjImg] = $subjImg"
echo "[bval_file] = $bval_file"
echo "[bvec_file] = $bvec_file"
echo "[outdir] = $outdir"
echo "[brainmask] = $brainmask"
echo "[round_bvals_1000] = $round_bvals_1000"
echo "[use_single_shell] = $use_single_shell"
echo "[del_files] = $del_files"

#####
# Implement QC
#####

## create intermediate file directory and QC summary text file
if [ ! -d $outdir/QC_output ]; then mkdir -p $outdir/QC_output;fi
subjQCDIR=$outdir/QC_output

echo "--- $subj QC summary ---" > $subjQCDIR/${subj}_QC.txt

## displacement (subject movement QC)
echo "--- Step 1: Calculating inter-volume motion ---"
mcflirt -in $subjDIR/$subjImg -refvol 0 -cost mutualinfo -mats -plots -rmsrel -rmsabs -o $subjQCDIR/${subjImg}_mcf
echo "mean relative motion: $(cat $subjQCDIR/${subjImg}_mcf_rel_mean.rms) " >> $subjQCDIR/${subj}_QC.txt
echo "mean absolute motion: $(cat $subjQCDIR/${subjImg}_mcf_abs_mean.rms) " >> $subjQCDIR/${subj}_QC.txt

## sum of squared differences (distortion QC)
# image intensity sum of squares (image intensity for each shell - mean intensity of each shell)^2
echo "--- Step 2: Calculating sum of squares for image intensity ---"
# find all unique b values
IFS=" " read -a bval_arr < "$subjDIR/$bval_file"

# round off bvals to nearest 1000 if needed
if [[ $round_bvals_1000 == 'y' ]]; then
    for idx in ${!bval_arr[@]}; do
        # rounded value = ((bval + 500) / 1000) * 1000
        bval_arr[$idx]=$(( $(( $(( ${bval_arr[$idx]} + 500 )) / 1000 )) * 1000 ))
    done
fi

# find number of shells
uniq_b=$(echo ${bval_arr[@]} | tr ' ' '\n' | sort -u | sort -n )
echo ${uniq_b[@]}

for search_b in ${uniq_b[@]}; do
    # exit if brain mask needs to be defined and no b0 values are there
    if [[ ${brainmask} == 'NIL' ]] && [[ ! $search_b == 0 ]]; then
        echo "ERROR: cannot generate brain mask, no b0 acquisitions detected"
    fi
    for idx in ${!bval_arr[@]}; do
        if [[ ${bval_arr[$idx]} == ${search_b} ]]; then
            fslroi $subjDIR/${subjImg}.nii $subjQCDIR/${subjImg}_b${search_b}_${idx} ${idx} 1
        fi;
    done
    fslmerge -t $subjQCDIR/${subjImg}_b${search_b} $subjQCDIR/${subjImg}_b${search_b}_*
    rm $subjQCDIR/${subjImg}_b${search_b}_*.nii.gz

    ## SNR (image quality QC)
    # mean signal from b0 across standard deviation
    # calculate SNR from b0 images
    if [[ ${search_b} == 0 ]]; then
    echo "--- Step 2.5: Calculating b0 SNR ---"
        if [[ ${brainmask} == 'NIL' ]]; then
            echo "Generating brain mask..."
            bet $subjQCDIR/${subjImg}_b0 $subjQCDIR/${subjImg}_bet_b0 -m
            brainmask=$subjQCDIR/${subjImg}_bet_b0_mask
        fi
        fslmaths $subjQCDIR/${subjImg}_b0 -Tmean -mas $brainmask $subjQCDIR/${subjImg}_mean_b0
        fslmaths $subjQCDIR/${subjImg}_b0 -Tstd -mas $brainmask $subjQCDIR/${subjImg}_std_b0
        fslmaths $subjQCDIR/${subjImg}_mean_b0 -div $subjQCDIR/${subjImg}_std_b0 $subjQCDIR/${subjImg}_SNR_b0
        echo "b0 data SNR: $(fslstats $subjQCDIR/${subjImg}_SNR_b0 -M)" >> $subjQCDIR/${subj}_QC.txt
    fi

    fslmaths $subjQCDIR/${subjImg}_b${search_b} -Tmean -mas $brainmask $subjQCDIR/${subjImg}_mean_b${search_b}
    fslmaths $subjQCDIR/${subjImg}_b${search_b} -sub $subjQCDIR/${subjImg}_mean_b${search_b} \
    -sqr $subjQCDIR/${subjImg}_SSD_b${search_b}
    echo "b${search_b} mean SSD: $(fslstats $subjQCDIR/${subjImg}_SSD_b${search_b} -M)" >> $subjQCDIR/${subj}_QC.txt
done

## residuals (tensor fit QC)
# use sse file
echo "--- Step 3: Calculating residuals ---"
if [ ! -d $outdir/fdt ]; then mkdir -p $outdir/fdt;fi
subjfdtDIR=$outdir/fdt
# run fsl dtifit
if [[ ! $use_single_shell == 'NIL' ]]; then
    # generate bvals and bvecs for 1000 shell
    echo "Extracting single shell for dtifit"
    # generate image with only b0 and desired shell
    count_frames=0
    for idx in ${!bval_arr[@]}; do
        if [[ ${bval_arr[$idx]} == 0 ]] || [[ ${bval_arr[$idx]} == $use_single_shell ]]; then
            echo "Image[$idx] = ${bval_arr[$idx]}"
            if [[ $count_frames == 0 ]]; then
                fslroi $subjDIR/${subjImg}.nii $outdir/fdt/${subj}_single_shell_${use_single_shell} ${idx} 1
                ((count_frames++))
            else
                fslroi $subjDIR/${subjImg}.nii $outdir/fdt/${subj}_single_shell_slice_${idx} ${idx} 1
                fslmerge \
                -t $outdir/fdt/${subj}_single_shell_${use_single_shell} $outdir/fdt/${subj}_single_shell_${use_single_shell} \
                $outdir/fdt/${subj}_single_shell_slice_${idx}
                ((count_frames++))
            fi
        fi
    done
    echo "Using $count_frames frames in total"
    rm $outdir/fdt/${subj}_single_shell_slice_*.nii.gz
    # matlab script should be in same directory as diffusion_QC.sh
    matlab -nodesktop -nosplash -nodisplay -r "addpath('$scriptdir/utilities'); \
    CBIG_DiffProc_extract_single_shell('$outdir/fdt','$subj', \
     '$subjDIR/$bval_file', '$subjDIR/$bvec_file', $use_single_shell, '$round_bvals_1000'); exit"
    dtifit -k $outdir/fdt/${subj}_single_shell_${use_single_shell} -o $subjfdtDIR/dti \
    -m $brainmask -r $outdir/fdt/${subj}_single_shell_${use_single_shell}.bvec \
    -b $outdir/fdt/${subj}_single_shell_${use_single_shell}.bval --sse
    echo "DTI tensor SSE: $(fslstats $subjfdtDIR/dti_sse -M)" >> $subjQCDIR/${subj}_QC.txt
else
    dtifit -k $subjDIR/$subjImg.nii -o $subjfdtDIR/dti -m $brainmask -r $subjDIR/$bvec_file -b $subjDIR/$bval_file --sse
    echo "DTI tensor SSE: $(fslstats $subjfdtDIR/dti_sse -M)" >> $subjQCDIR/${subj}_QC.txt
fi

## remove QC files and only leave txt file(to save space)
if [[ $del_files == 'y' ]]; then
    #rm $subjDIR/${subjImg}.nii
    if [[ -d $subjQCDIR ]]; then
        cd $subjQCDIR
        ls -I *.txt | xargs rm -r
    fi
fi

