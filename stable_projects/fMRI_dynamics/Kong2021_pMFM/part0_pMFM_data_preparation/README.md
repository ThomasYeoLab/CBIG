# Parcellated TC and FCD generation for both Desikan and Schaefer 100 parcellations
* `PARCELLATION_NAME/input` this folder contains the data to perform the estimation process of pMFM. Files in this folder are used by the scripts in `scripts` folder.
* `Desikan/scripts` this folder contains the scripts to generate TC of fMRI signal in Desikan parcellation
    * `CBIG_pMFM_step1_generate_TC_desikan.m` this file is the wrapper to generate the TC in Desikan parcellation.
    * `CBIG_pMFM_step2_generate_FCD_desikan.m` this file is the wrapper to generate the FCD matrices.
    
* `Schaefer100/scripts` this folder contains the scripts to generate TC of fMRI signal in Schaefer 100 parcellation
    * `CBIG_pMFM_step1_generate_TC_desikan.m` this file is the wrapper to generate the TC in Schaefer 100 parcellation.
    * `CBIG_pMFM_step2_generate_FCD_desikan.m` this file is the wrapper to generate the FCD matrices.
        
    
# Usage
* Go to working directory using command: `cd $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Kong2021_pMFM/part0_pMFM_data_preparation/PARCELLATION_NAME/scripts` and follow the step order from step 1 to step 3 to run the scripts.
* For Matlab scripts, you can run it directly in terminal.