## Background

This folder contains the modified TBSS scripts to run using CSCHPC's cluster. For full details on how to run TBSS, please visit FSL's user 
guide for using TBSS [https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TBSS/UserGuide](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TBSS/UserGuide).

### Usage

TBSS comprises of 4 main steps. 

- Preparation prior to TBSS

  Prior to TBSS, each subject's FA volume needs to be put into a single directory. Any secondary metrics that you want to project to the TBSS skeleton
  has to be put in a folder with the metric name in the main directory. The name of the secondary metric image must be the same as the FA image. Here is
  an example where we wish to project the FA images and MD images of 3 subjects to the TBSS skeleton:

```
TBSS output directory 
|-- sub-1_FA.nii.gz 
|-- sub-2_FA.nii.gz 
|-- sub-3_FA.nii.gz 
|-- MD 
    |-- sub-1_FA.nii.gz
    |-- sub-2_FA.nii.gz 
    |-- sub-3_FA.nii.gz
```



- Step 1: Preprocessing

  Each image is moved to a FA folder and a mask is created for each subject. Start this step by changing directory to the folder with the FA volumes 
  and calling the following command:
```
export scriptdir=$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_diffusion_processing2022/TBSS
export tbssdir=output_dir
$script_dir/tbss_1_preproc_torque *.nii.gz
```

- Step 2: Registration

  Each image is registered to a template.
```
export scriptdir=$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_diffusion_processing2022/TBSS
export tbssdir=output_dir
$script_dir/tbss_2_reg_torque -T
```

- Step 3: Post-registration
  
  The TBSS skeleton is formed from the FA images.
```
export scriptdir=$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_diffusion_processing2022/TBSS
export tbssdir=output_dir
$script_dir/tbss_3_postreg_torque -S
```

- Step 4: Post-registration

  A final thresholding is applied to the TBSS skeleton.
```
tbss_4_prestats 0.2
```

- Step 5: Projection of other metrics (Optional)

  Other metrics can be projected to the TBSS skeleton. Make sure that they are saved in the main directory with the format specified above.
```
export scriptdir=$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_diffusion_processing2022/TBSS
export tbssdir=output_dir
$script_dir/tbss_non_FA_torque MD
```

----

## Bugs and Questions

Please contact Leon Ooi at leonooiqr@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.