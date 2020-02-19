# Data Release
In this folder, we release the factor compositions of the ASD participants included in our paper. We also release the factor-specific hypo/hyper RSFC patterns for the three-factor model.

----

## Reference

Siyi Tang*, Nanbo Sun*, Dorothea L. Floris, Xiuming Zhang, Adriana Di Martino, B.T. Thomas Yeo. [Reconciling Dimensional and Categorical Models of Autism Heterogeneity: a Brain Connectomics & Behavioral Study](https://doi.org/10.1016/j.biopsych.2019.11.009). Biological psychiatry, in press.

----

## Looking For The Factor-Specific RSFC Patterns?

Factor-specific hypo/hyper RSFC patterns (i.e.,E(RSFC patterns|Factor)) for the three-factor model are located at `files_polarLDA`.

----

## Want To Infer Factor Compositions of New Participants?
If you have RSFC data (419 x 419 matrices) of new participants and want to infer the factor compositions of these participants with our estimated model parameters, we have provided the following files for you:
1. `files_FC2doc/ref_mean_std.mat` contains the mean and standard deviation of RSFC data of control participants in ABIDE-II+GENDAAR after regressing out nuisance variables. You need to perform z-normalization of your new participant's RSFC data with respect to this mean & standard deviation before inferring the factor compositions.
2. `files_polarLDA` folder contains the estimated model parameters: `final.beta`, `final.rho`, `final.gamma` & `final.other` for the three-factor model. These parameters were estimated from ASD participants in ABIDE-II+GENDAAR, and will be used to infer factor compositions of new participants.

For more details on how to call the relevant functions to infer factor compositions of new participants, please refer to the example `../examples`.

----

## Factor Compositions of ASD Participants In Our Study


Factor compositions (three-factor model) of the ASD participants included in our paper are released as a spreadsheet in `factorCompositions_ASD.csv`. Specifically, the spreadsheet includes factor compositions (i.e., Pr(Factor|Participant)) of the 306 ASD participants in ABIDE-II+GENDAAR datasets as well as 166 ASD participants in ABIDE-I dataset.

#### Column Definitions
The spreadsheet `factorCompositions_ASD.csv` has 6 columns: `SITE_ID`, `REPOSITORY`, `SUB_ID`, `FACTOR1_PROB`, `FACTOR2_PROB` & `FACTOR3_PROB`.
1. `SITE_ID`: This column corresponds to the acquisition sites of the participants.
2. `REPOSITORY`: This column corresponds to the data repository that the participants belong to, i.e., ABIDE-I, ABIDE-II or GENDAAR.
3. `SUB_ID`: This column corresponds to the subjects' IDs defined in the respective data repository.
4. `FACTOR1_PROB`: This column corresponds to the probability of factor 1 (i.e., Pr(Factor 1|Participant)).
5. `FACTOR2_PROB`: This column corresponds to the probability of factor 2 (i.e., Pr(Factor 2|Participant)).
6. `FACTOR3_PROB`: This column corresponds to the probability of factor 3 (i.e., Pr(Factor 3|Participant)).
