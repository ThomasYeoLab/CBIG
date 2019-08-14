This folder contains scripts about how to replicate all the findings from our Biological Psychiatry paper:

(1) PLS analysis 
(2) Cross-validation
(3) Validation on task FC data

Note that all filenames and directories in these scripts **only work for the CBIG lab**.

----

## Data
The scripts uses data from the UCLA Consortium for Neuropsychiatric Phenomics (CNP). Preprocessed data can be found in 

```
/mnt/eql/yeo9/data/UCLAconsortiumNeuropsych
```

The data can be (re-) downloaded at https://openneuro.org/datasets/ds000030/versions/00016

----

## Code

The wrapper function is `CBIG_VK2019_replication_wrapper.m`.

The main function (`PLS_analysis.m`) is located in

``` 
$CBIG_CODE_DIR/external_packages/matlab/non_default_packages/PLS_MIPlab
```


and was modified from the code written by Dimitri Van De Ville, Daniella Zoller and Valeria Kebets, available at https://miplab.epfl.ch/index.php/software/PLS

Results can be compared with `PLSresults.mat`, `PLSresults_1000permuts.mat`, `PLS_bootstrapLoadings_500bootstraps.mat`, `PLS_bootstrapResults_500bootstraps.mat`, `PLS_5fold_crossval_1000permuts.mat`, and `PLS_taskValidation_1000permuts.mat` in 

```
/mnt/eql/yeo1/CBIG_private_data/replication/stable_projects/disorder_subtypes/Kebets2019_TransdiagnosticComponents/results
``` 

----
