## Reference

Jianxiao Wu, Gia H. Ngo, Alexander Schaefer, Douglas Greve, Jingwei Li, Tong He, Bruce Fischl, Simon B. Eickhoff, B.T. Thomas Yeo. [**Accurate Nonlinear Mapping between MNI152/Colin27 Volumetric and FreeSurfer Surface Coordinate Systems**](http://people.csail.mit.edu/ythomas/publications/2018VolSurfMapping-HBM.pdf), *Human Brain Mapping*, 2018.

----

## Background

The results of most neuroimaging studies are reported in volumetric (e.g., MNI152) or surface (e.g., fsaverage) coordinate systems. Accurate mappings between volumetric and surface coordinate systems can facilitate many applications, such as projecting fMRI group analyses from MNI152/Colin27 to fsaverage for visualization, or projecting resting-state fMRI parcellations from fsaverage to MNI152/Colin27 for volumetric analysis of new data. 

In this project, we evaluated three approaches for mapping data between MNI152/Colin27 and fsaverage coordinate systems. Two of the approaches (MNIsurf/Colinsurf and Affine) are currently widely used. A third approach (registration fusion) was previously proposed, but not widely adopted. Two implementations of the registration fusion (RF) approach were considered, with one implementation utilizing the Advanced Normalization Tools (ANTs). We found that RF-ANTs performed the best in general. 

The visualisation of three anatomical regions (represented as probabilistic maps) projected from MNI152 to fsaverage is shown below. The black boundaries show the "ground truth" parcellation. 

<img src="bin/images/root_readme_img.png" height="500" />

----

## Stand-alone Usage

To use volume-to-fsaverage mappings without the trouble of downloading our entire repository, just download the `bin/final_warps_FS5.3` and `bin/scripts_stand_alone_for_MNI_fsaverage_projection` (or `bin/scripts_stand_alone_for_Colin_fsaverage_projection`) folders. Put them in the same directory; then follow the README.md in `bin/scripts_stand_alone_for_MNI_fsaverage_projection` (or `bin/scripts_stand_alone_for_Colin_fsaverage_projection`).

Note that FreeSurfer and Matlab need to be installed before the stand-alone scripts can be run.

----

## Data Release

The folder `bin` contains data files necessary for this project, as well as the final mappings released.

- `bin/GSP_subjectid.csv` is a file listing all 1490 GSP subject IDs used in this project

- `bin/liberal_cortex_masks_FS5.3` folder: contains the liberal cortical masks of MNI152 and Colin27, which are used in fsaverage-to-volume atlas projections

- `bin/final_warps_FS5.3` folder: contains the MNI152-fsaverage and Colin27-fsaverage mappings generated in this project, for both RF-M3Z and RF-ANTs approaches, using 1490 GSP subjects. These mappings can be used by calling scripts from `bin/scripts_final_proj_FS5.3`

----

## Code Release

- `registration_fusion` folder: contains codes for running RF approaches (RF-M3Z and RF-ANTs) from scratch. Specifically, `registration_fusion/scripts_vol2surf` contains codes for obtaining volume-atlas-to-fsaverage mapping, while `registration_fusion/scripts_surf2vol` contains codes for obtaining fsaverage-to-volume-atlas mapping. See README.md in the respective folders for instructions on implementation.

- `bin/scripts_final_proj_FS5.3` folder: contains codes for  using final warps from RF-approaches (RF-M3Z and RF-ANTs) directly. Refer to README.md in the folder for instructions on implementation.

- `freesurfer_baseline` folder:  contains codes for running Affine or MNIsurf/Colinsurf approach. Refer to README.md in the folder for instruction on implementation.

Note that this project uses generic functions from other folders, which may be updated over time. To download the version of the code that was last tested, you can either

- visit this link: https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.5.0-Wu2017_RegistrationFusion

or

- run the following command, if you have Git installed
```
git checkout -b Wu2017_RegistrationFusion v0.5.0-Wu2017_RegistrationFusion
```

----

## Updates

- Release v0.5.0 (15/04/2018): Initial release of Wu2017 Registration Fusion project.
- Release v0.6.4 (13/08/2018): Added stand-alone scripts for Colin27-to-fsaverage projections; added example.

----

## Bugs and Questions

Please contact Jianxiao Wu at vesaveronica@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
