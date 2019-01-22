This README includes instruction for unit-testing a coordinate-based meta-analysis with the author-topic model.

## Reference
Gia H. Ngo, Simon B. Eickhoff, Minh Nguyen, Peter T. Fox,  R. Nathan Spreng, B. T. Thomas Yeo. [Beyond Consensus: Embracing Heterogeneity in Neuroimaging Meta-Analysis](https://www.biorxiv.org/content/early/2017/06/19/149567). BioRxiv preprint

----

## Data
For this unit test, we used activation coordinates from the self-generated thought dataset included under `../SelfGeneratedThought/MNI152_ActivationCoordinates`.

----

## Run

In Matlab, run `CBIG_AuthorTopic_UnitTest.m`

If the local environment has been set up properly (see the project's main [README](https://github.com/ThomasYeoLab/CBIG/stable_projects/meta-analysis/Ngo2019_AuthorTopic)), the unit test would perform the following steps:

1. Prepare input data as a Matlab `.mat` file in the format required by the inference functions and save it under `./data` directory.
2. Perform inference to estimate the model parameters. For testing, inference will only be performed up to 3 cognitive components with 3 re-initializations each. All output files are saved under `./outputs` directory.
3. Visualize components of the 2-component solution. The output figures will be saved under `./figures` directory.
4. Compute Bayesian Information Criterion (BIC) across different number of components. The output files and figures will be saved under `./BIC` directory.
