Diffusion processing project depends on real diffusion MRI data, which may vary from case to case. Please follow the instructions below to run each example:

## diffusionQC
Please refer to the top level `CBIG2022_DiffProc` folder for specific instructions relating to this script.
### Preparing the example data
Data is from group 0 of the NKI diffusion data, which can be accessed [here](http://fcon_1000.projects.nitrc.org/indi/enhanced/data/download_nii8.html).
We use participant `A00059845` for this example. 
To prepare the data, create an `input` folder and place each participant's data inside (e.g. `input\A00059845`). The diffusion image, bval and bvec files should be in the same folder.
### Running the script
To run the script, the appropriate arguments have to be provided to `CBIG_DiffProc_diffusionQC.sh`. 
You can use `CBIG_DiffProc_diffusionQC_run_example.m` as a guide on how to call the script. 
You may check the results of your processing is consistent with the example using `CBIG_DiffProc_diffusionQC_check_example_results.m`.

## TBSS
Please refer to the `TBSS` folder for specific instructions relating to this script.
### Preparing the data
Data is from group 0 and group 1 of the NKI diffusion data, which can be accessed [here](http://fcon_1000.projects.nitrc.org/indi/enhanced/data/download_nii8.html).
We use participant `A00059845` and `A00058145` for this example. Please generate the diffusion-tensor data for these 2 participants using the diffusionQC pipeline from above.
Place the FA images from the `fdt` folder from the diffusionQC pipeline into an output folder. In the same folder create a sub-folder named `MD` and place the MD images from the same `fdt` folder.
Rename the MD images in the `MD` to match the names of the FA images in the output folder (Might be confusing but it is necessary to run the pipeline!). See the folder structure below as an example.
```
TBSS output directory 
|-- A00059845_FA.nii.gz 
|-- A00058145_FA.nii.gz 
|-- MD 
    |-- A00059845_FA.nii.gz
    |-- A00058145_FA.nii.gz 
```
### Running the script
To run the script, ensure you have FSL installed, and call them in the right order - refer to their user guide [here](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TBSS/UserGuide). 
You can use `CBIG_DiffProc_TBSS_run_example.m` as a guide on how to call the script. 
You may check the results of your processing is consistent with the example using `CBIG_DiffProc_TBSS_check_example_results.m`.
Additional Note: The script has been run without changing the default thresholds - the TBSS skeleton could be further improved by tuning the registration and the FA thresholds.

## AMICO
Please refer to the `AMICO` folder for specific instructions relating to this script.
### Preparing the data
Data is from the NODDI example dataset found [here](https://github.com/daducci/AMICO/wiki/NODDI). 
To run the script successfully, the file type needs to be changed to nifti, and the file names will have to be changed to match the script.
Please follow the suggested data structure outlined in the `AMICO` folder. To change the file type from the img and hdr format to nifti, you can use the command `fslchfiletype NIFTI_GZ name_of_file`.
For the example data in particular, you can use `fslchfiletype NIFTI_GZ subject_01`.
The `AMICO_example_data` folder has a template of the folder structure that you may follow. The data inside the folders have been removed to conserve space on github.
### Running the script
You can use `CBIG_DiffProc_AMICO_run_example.m` as a guide on how to call the script. 
You may check the results of your processing is consistent with the example using `CBIG_DiffProc_AMICO_check_example_results.m`.

## MRtrix
There is no example reference for MRtrix processing as the resulting files are too big to be stored on Github. The example reference is stored on the CBIG HPC servers and can only be run within the CBIG environment.