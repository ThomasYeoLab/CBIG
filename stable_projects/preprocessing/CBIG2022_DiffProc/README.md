## Reference

Please visit the various subfolders for the references for each pipeline.

----

## Background

This folder contains a diffusion processing pipeline written by CBIG group. The pipelines include:
- Diffusion quality checks (QC)
- Tract based spatial statistics (TBSS) skeleton generation
- AMICO (for generating NODDI images)
- MRtrix tractography

----

## Code Release

### Download stand-alone repository

Since the whole GitHub repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link:
[https://github.com/ThomasYeoLab/Standalone_CBIG2022_DiffProc](https://github.com/ThomasYeoLab/Standalone_CBIG2022_DiffProc)


### Download whole repository

Except for this project, if you want to use the code for other stable projects from our lab as well, you need to download the whole repository.

- To download the version of the code that was last tested, you can either

  - visit this link:
  [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.29.8-CBIG2022_DiffProc-updates](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.29.8-CBIG2022_DiffProc-updates)
  
  or
  
  - run the following command, if you have Git installed
  
  ```
  git checkout -b CBIG2022_DiffProc-updates v0.29.8-CBIG2022_DiffProc-updates
  ```

### Usage 

- How to run our processing pipeline?

  The current processing scripts only work with minmally processed data (motion correction, dark slice replacement, distortion correction and eddy current correction will have to be completed beforehand). 
  We recommend running `CBIG_diffusionQC.sh` to check the quality of your diffusion data as a first step. The script fits a DTI model to the data for each subject, the resulting images can be used for 
  secondary processing steps such as TBSS. 

- Quality checks
  This script generates the quality assessment metrics in a text file. The subject ID, location and names the diffusion files will have to be provided. An example of calling the script is provided below:
```
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG2022_DiffProc/CBIG_diffusionQC.sh \
    --subj $subject --subjDIR $projdir/$subject/T1w/Diffusion \
    --subjImg data.nii.gz --outdir $base_outdir/S1200/individuals/$subject \
    --bval_file bvals --bvec_file bvecs \
    --brainmask $projdir/$subject/T1w/Diffusion/nodif_brain_mask.nii.gz 
```
  In this function, the following metrics are calculated:
  1. Intervolume motion
  2. b=0 images SNR
  3. Sum of squared differences of signal intensity for each diffusion shell
  4. Sum of residuals over all voxels when fitting to the DTI model
  
  These results are saved in a text file which the user may use to impose a selection criteria to exclude subjects whose diffusion images are of poor quality.
  The limits should differ from study to study, as the acquisition protocol can vary widely. However, as an example, a user may choose to exclude the subject 
  if he/she exceeds any of the following criteria `motion > 1mm`, `b0 SNR < 20`, `Sum of squared differences > 2e5`, `Sum of residuals > 22`. 

- TBSS skeleton generation

  A folder containing all fractional anisotropy (FA) volumes of the subjects to process TBSS will have to be provided. For more details about this pipeline, refer to the readme in the TBSS folder.

- AMICO pipeline for multishell diffusion data

  A folder containing the diffusion images of each subject will have to be provded. For more details about this pipeline, refer to the readme in the AMICO folder.

- MRtrix structural connectivity matrix generation

  Each subject needs to have a recon-all folder, diffusion folder. Parcellations of each subject will need to be converted to the native volume space. For more details about this pipeline, refer to the readme in the MRtrix folder.

----

## Updates

- Release v0.29.8 (6/3/2024): Update to unit tests and examples.
- Release v0.29.6 (16/1/2024): Update to include examples from publicly available data.
- Release v0.25.0 (23/9/2022): Initial release of CBIG diffusion processing pipeline.
  
----

## Bugs and Questions

Please contact Leon Ooi at leonooiqr@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
