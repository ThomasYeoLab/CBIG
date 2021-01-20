# Reference

Kebets V, Holmes AJ, Orban C, Tang S, Li J, Sun N, Kong R, Poldrack RA, Yeo BTT. **Somatosensory-Motor Dysconnectivity Spans Multiple Transdiagnostic Dimensions of Psychopathology**. Biological Psychiatry, 86(10): 779-791, 2019

----

# Background

The societal and economic burden of mental illness is increasingly recognized, as 1 in 7 people in Singapore has experienced a psychiatric disorder in their lifetime. The biology underlying psychiatric disorders remain poorly understood. Two patients with different psychiatric diagnoses often present the same symptoms (e.g., depressed mood) and cognitive deficits (e.g., memory problems), suggesting that psychiatric disorders are not distinct categories, but instead share common dysfunctional neurobiology. We found shared patterns of brain dysfunction across healthy individuals and individuals with bipolar disorder, schizophrenia and attention-deficit/hyperactivity disorder (ADHD). The magnitude of brain dysfunction was associated with symptom severity and cognitive impairment. Moreover, brain regions involved in motor control were particularly disturbed, and might be related to the motor deficits that are commonly experienced by patients with psychiatric disorders. These findings give new insights on potential mechanisms underlying symptoms and cognitive deficits that are shared among psychiatric disorders.

----

# Important Note
- We no longer support this project.
- The last supported git repo version is v0.17.0-Fix_Absolute_Path. To access this:
    - visit this link:
    [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.17.0-Fix_Absolute_Path](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.17.0-Fix_Absolute_Path)

    or

    - run the following command, if you have Git installed
 
    ```
    git checkout -b Kebets2019_TransdiagnosticComponents v0.17.0-Fix_Absolute_Path
    ```
- Last working config file is the CBIG_sample_config.sh in the last supported git repo version.

----

# Data Release
The main findings of this study are available in the `data_release` folder. Specifically, the data include:
1. RSFC data of all subjects included in the analysis
2. RSFC and behavioral composite scores of all subjects associated with the 3 transdiagnostic components
3. RSFC and behavioral loadings associated with each component
4. Supplementary Table S7 (complete table)

The subjects included in the analysis (N=224) and their primary diagnosis are listed in `LC_composite_scores.csv`.

----

# Code Release
### Download stand-alone repository
Since the whole Github repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link: [https://github.com/ThomasYeoLab/Standalone_Kebets2019_TransdiagnosticComponents](https://github.com/ThomasYeoLab/Standalone_Kebets2019_TransdiagnosticComponents)

### Download whole repository
Except for this project, if you want to use the code for other stable projects from out lab as well, you need to download the whole repository.

- To download the version of the code that was last tested, you can either

    - visit this link:
    [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.17.0-Fix_Absolute_Path](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.17.0-Fix_Absolute_Path)

    or

    - run the following command, if you have Git installed
 
    ```
    git checkout -b Kebets2019_TransdiagnosticComponents v0.17.0-Fix_Absolute_Path
    ```

----

# Usage
### Set up
1. Make sure you have installed: Matlab 2018b, Freesurfer 5.3.0, Python 3.6
2. Follow `$CBIG_CODE_DIR/setup/README.md` to setup the config file. Instead of using `CBIG/setup/CBIG_sample_config.sh`, 
you need to use `$CBIG_CODE_DIR/stable_projects/disorder_subtypes/Kebets2019_TransdiagnosticComponents/replication/config/CBIG_VK2019_tested_config.sh`.
3. To create a Python environment similar to what was used for the project, install Anaconda with Python 3.6, and create the Anaconda environment from `config/CBIG_VK2019_python_env.yml` file by running `conda env create -f config/CBIG_VK2019_python_env.yml`

### Unit tests
The unit test is mainly for internal use.

### Examples
In the folder `examples`, we provide a toy example code on how to derive & visualize latent components. The user can run this example to verify the set up is correct. 

### Replication
In the folder `replication`, we provide scripts to replicate the main findings of our paper.
The code for the PLS analysis is located in

``` 
$CBIG_CODE_DIR/external_packages/matlab/non_default_packages/PLS_MIPlab
```

and was modified from the code written by Valeria Kebets, Daniela Zöller and Dimitri Van De Ville (https://github.com/danizoeller/myPLS).


----

# Updates

- Release v0.11.0 (14/06/2019): Initial release of the Kebets2019_TransdiagnosticComponents project

- Release v0.13.3 (18/07/2019): Release of code for Kebets2019_TransdiagnosticComponents

- Release v0.13.4 (14/08/2019): Updated `Table S7_sumCorr.csv` & `Schaefer2018_400Parcels_17Networks_19Subcortical.csv` with corrected labels

- Release v0.14.2 (09/09/2019): Corrected group resampling in bootstrap function, and added random number generation to bootstrap & permutation functions

- Release v0.15.1 (08/10/2019): Released RSFC matrices, and updated `Table S7_sumCorr.csv` & `Schaefer2018_400Parcels_17Networks_19Subcortical.csv` with corrected labels

- Release v0.17.0 (19/02/2020): Added new environment variables to avoid possible problems caused by hard-coded absolute paths, added names of behavioral variables, and reduced the size of the example data.
    
- Release v.x.x (15/05/2020): Split the example wrapper into functions to generate and check the example data.


----

# Bugs and Questions
Please contact Valeria Kebets at valkebets@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

