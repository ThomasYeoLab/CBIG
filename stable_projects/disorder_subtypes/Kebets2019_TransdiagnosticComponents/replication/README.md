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

The main function (i.e., PLS analysis) was modified from the code written by Dimitri Van De Ville, Daniella Zoller and Valeria Kebets, available at [https://miplab.epfl.ch/index.php/software/PLS].

Results can be compared with `PLSresults_all54Scores_N224.mat` in 

```
/mnt/eql/p1/users/external/vkebets/CBIG_private/stable_projects/disorder_subtypes/Kebets2019_TransdiagnosticComponents/replication/correct_output
```

----
