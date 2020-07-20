# Spatial Distortion Correction (SDC) README

## Overview

**Spatial distortion correction (SDC) is one of the preprocessing steps that aim to remove the distortion in EPI functional images caused by magnetic field inhomogeneities.** These magnetic field inhomogeneities are due to different magnetic susceptibilities of tissues within the head. The results of these magnetic field inhomogenieties include spatial distortions and signal loss, particularly in inferior frontal and temporal regions. SDC utilizes the information of **fieldmaps** to unwarp the distorted functional images. If you consider running SDC as one of the preprocessing steps, please note that SDC should be run after motion correction and before registration step (See `example_with_distortion_correction_config.txt`for an example configuration file).

## Preparation

In order to run SDC, please make sure that the following images and information are available:

* EPI functional image
* Fieldmaps (as described in the next section)
* Echo time (TE) of the functional image (in ms)
* Effective echo spacing (EES) of the fucntional image (in ms)
* Motion correction transformation matrices

TE and EES of the functional image can be typically found in the JSON file of the associated functional image if the dataset follows BIDS format.

Motion correction transformation matrices should be available if motion correction step is run successfully. For our CBIG preprocessing pipeline, the transformation matrices file is named as `<subject_id>/bold/<run_number>/<bold_stem>_mc.cat`.

## Fieldmaps

There are two main types of fieldmaps, the fieldmaps are typically located under `/fmap` directoies if the dataset follows Brain Imaging Data Structure (BIDS) format. These two types of fieldmaps are as following:

* Two magnitude images with different echo times and one phase-difference image (**mag+phasediff**).
* Two images with opposite phase encoding directions (**oppo_PED**). Note that we only support 'j+' and 'j-' phase encoding directions for now.

## How to run the SDC script?

### mag+phasediff

If the fieldmaps include 2 magnitude images and 1 phase difference image, besides the above-mentioned images and information (as described in Preparation section), you will also need the difference in echo time between 2 magnitude images (in ms). As for the magnitude image, you can choose either one as the input image. An example command is shown as the following:

`./CBIG_preproc_spatial_distortion_correction -s sub-NDARAA536PTU -d ~/storage/fMRI_preprocess -bld '001' -BOLD_stem _rest -fpm "mag+phasediff" -m /data/HBN/rawData_release1_4/SI/sub-NDARAA536PTU/fmap/ sub-NDARAA536PTU_magnitude1.nii.gz -p /data/HBN/rawData_release_1_4/SI/sub-NDARAA536PTU/fmap/sub-NDARAA536PTU_phasediff.nii.gz -delta 4.76 -ees 0.55 -te 40`

More specifically,

* `-fpm` is for fieldmap processing method, the processing method here is "mag+phasediff"
* `-m` is the path to one of the magnitude images
* `-p` is the path to the phase difference image
* `-delta` is the absolute difference in echo time between 2 magnitude images (in ms)
* `-ees` is the effective echo spacing of the functional image (in ms)
* `-te` is the echo time of the functional image (in ms)

### oppo_PED

If the fieldmaps are in opposite phase encoding directions,  besides the above-mentioned images and information (as described in Preparation section), you will also need the total readout time (TRT) of 2 fieldmaps (in seconds). 

**[IMPORTANT]**: Please note that we only support 'j+' and 'j-' phase encoding directions for now. To find out what are the phase encoding directions of the fieldmaps, please check the orientation of the fieldmaps by running `mri_info <fielmap>`. For example, if the orientation is RAS, then it means Anterior (also Right and Superior) is the postive directions. Therefore, under RAS orientation, if a fieldmap has a phase encoding direction of PA, then it corresponds to 'j+' as the voxel coordinate increases from Posterior to Anterior. On the other hand, if a fieldmap has a phase encoding direction of AP, then it corresponds to 'j-' as the voxel coordinate decreases from Anterior to Posterior. Double check with the JSON files associated with the fieldmaps if the dataset follows BIDS format.

An example command is shown as the following:

`./CBIG_preproc_spatial_distortion_correction.csh -s sub-NDARWN691CG7 -d ~/storage/fMIR_preprocess -bld '001 002' -BOLD_stem _rest -fpm "oppo_PED" -j_minus /data/HBN/rawData_release1_4/RU/sub-NDARWN691CG7/fmap/sub-NDARWN691CG7_dir-AP_acq-fMRI_epi.nii.gz -j_plus /data/HBN/rawData_release1_4/RU/sub-NDARWN691CG7/fmap/sub-NDARWN691CG7_dir-PA_acq-fMRI_epi.nii.gz -j_minus_trt 0.04565 -j_plus_trt 0.04565 -ees .580013000 -te 30.00`

More specifically,

* `-fpm` is for fieldmap processing method, the processing method here is "oppo_PED"
* `-j_minus` is the path to `j-` fieldmap
* `-j_plus` is the path to `j+` fieldmap
* `-j_minus_trt` is the total readout time of `j-` fieldmap (in seconds)
* `-j_plus_trt` is the total readout time of `j+` fieldmap (in seconds)
* `-ees` is the effective echo spacing of the functional image (in ms)
* `-te` is the echo time of the functional image (in ms)

Total readout time can be typically found in the JSON files of the associated fieldmaps if the dataset follows BIDS format.

## Quality Control (QC)

There are two ways to check if SDC has done its job. 

The qulitative way is to visually compare the functional images before and after SDC. The result after distortion correction should be more similar to the structrual image, especially in frontal and temporal regoins. 

The quantitative way, which is also the way we adopt for QC is to compare the BBR costs with and without SDC. By right, the BBR cost should be lower with SDC (See `quality_control_readme.md` for more details).
