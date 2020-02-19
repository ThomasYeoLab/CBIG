# Data Release
In this folder, we have made available the subjects' RSFC data, as well as their RSFC and behavioral composite scores for the 3 significant latent components (LCs 1-3) obtained with Partial Least Squares (PLS) in healthy individuals and individuals with ADHD, bipolar disorder, schizophrenia and schizoaffective disorder.
We also release the RSFC and behavioral loadings associated with LCs 1-3.

----

## Reference
Kebets V, Holmes AJ, Orban C, Tang S, Li J, Sun N, Kong R, Poldrack RA, Yeo BTT. **Somatosensory-Motor Dysconnectivity Spans Multiple Transdiagnostic Dimensions of Psychopathology**. Biological Psychiatry, 86(10): 779-791, 2019

----

## Subjects' RSFC data
The `RSFC_matrices` folder contains the RSFC data of all subjects included in the PLS analysis. 
The RSFC data are 419 x 419 matrices, with the first 400 ROIs being Schaefer's 400 cortical regions (Schaefer et al., 2018), and the last 19 ROIs being subcortical regions from the FreeSurfer segmentation (Fischl et al., 2002).
The spreadsheet `Schaefer2018_400Parcels_17Networks_19Subcortical.csv` lists the 419 ROIs (as ordered in the .mat files). 
The RSFC matrices can be plotted using

```
$CBIG_CODE_DIR/utilities/matlab/figure_utilities/CBIG_Plot_Schaefer400_17Networks19SubcorRearrCorrMat_WhiteGrid.m
``` 

Note that the ROIs are re-arranged for visualisation.

----

## PLS loadings
### RSFC and behavioral loadings of Latent Components 1-3
The RSFC/behavioral loadings (i.e., Pearson's correlations between the original RSFC/behavioral data and RSFC/behavioral composite scores) associated with LCs 1-3 can be found in the `PLS_loadings` folder.
The behavioral variables (ordered as the behavioral loadings) are listed in `behav_vars.csv`.

### Subjects' RSFC and behavioral composite scores of Latent Components 1-3
Subject-specific RSFC and behavioral composite scores for LCs 1-3 of the 224 participants included in our paper are released in the spreadsheet `/PLS_loadings/LC_composite_scores.csv`. 

##### Column definitions
The spreadsheet `LC_composite_scores.csv` has 8 columns: 
1. `SUB_ID`: Subjects' ID
2. `SUB_DIAGNOSIS`: Subjects' primary diagnosis, i.e., CONTROL, ADHD, BIPOLAR, SCHZ (schizophrenia) or SCHIZAFF (schizoaffective disorder)
3. `LC1_RSFC_COMPOSITE_SCORE`: Subjects' RSFC composite scores on the first latent component (LC1)
4. `LC2_RSFC_COMPOSITE_SCORE`: Subjects' RSFC composite scores on the second latent component (LC2)
5. `LC3_RSFC_COMPOSITE_SCORE`: Subjects' RSFC composite scores on the third latent component (LC3)
6. `LC1_BEHAV_COMPOSITE_SCORE`: Subjects' behavioral composite scores on the first latent component (LC1)
7. `LC2_BEHAV_COMPOSITE_SCORE`: Subjects' behavioral composite scores on the second latent component (LC2)
8. `LC3_BEHAV_COMPOSITE_SCORE`: Subjects' behavioral composite scores on the third latent component (LC3)
