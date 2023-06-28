
## Content

* Sub-folder `data/` contains all necessary data.

    ```
    1000subjects_clusters007_ref.mat
    Desikan34ROIsperHemisphere_regionDescription.mat
    Example_Estimated_Parameter                            
    ThomasYeo_52_components_read_out.mat                 
    ```
    
* Sub-folder `lib/` contains all functions.

   ```
   CBIG_MFMem_rfMRI_mfm_Num2Color.m
   CBIG_MEMem_rfMRI_FDR.m
   ```
   
*  `CBIG_MFMem_rfMRI_DrawBarchart_I_Desikan68_into_Yeo7_fsaverage5.m` 
   This function is to show mean value of model parameter I in Yeo's 7-network.
    
*  `CBIG_MFMem_rfMRI_DrawBarchart_W_Desikan68_into_Yeo7_fsaverage5.m` 
   This function is to show mean value of model parameter W in Yeo's 7-network.

*  `CBIG_MFMem_rfMRI_DrawOverlap_WI_Desikan68_Yeo7_fsaverage5.m` 
   This function is to show estimated W and I in Desikan parcellation with Yeo's 7-network boundary

*  `CBIG_MFMem_rfMRI_plot_corr_of_VonEconomo.m` 
   This function is to show correlation of estimated W and I in Desikan parcellation with VonEconomo cytoarchitecture data

----

## Usage

The folder contains all data, functions and scripts to plot the brain map and bar chart of the estimated parameters obtained in Step 1: Estimation.

```
CBIG_ReadNCAvgMesh.m
CBIG_DrawSurfaceMapsWithBoundary.m
```
Make sure that you are under `Wang2018_MFMem/step3_plot` to run the functions above.

----

## Bugs and Questions

Please contact Dr. Peng Wang (pengwanghome@gmail.com) and Prof. Dr. Thomas Yeo (yeoyeo02@gmail.com).
