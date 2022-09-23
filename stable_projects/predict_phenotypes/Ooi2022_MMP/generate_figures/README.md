# Regression codes
This folder contains code and text files required to generate the figures from the paper.

### Prerequisites
The code imports libraries from scipy, pandas, numpy and seaborn. If you do not have this in your python environment, you can use `../replication/config/Ooi2022_MMP_figures.yml` to set up your environment. 
Additionally, before saving the code, please save all regression results from kernel ridge regression, linear ridge regression and elastic net using the scripts `CBIG_MMP_<dataset>_collate_result_wrapper.m`. This would have resulted in a .mat file saving the results of each single-feature-type model result for each component score and the "grand average" (Average over all original scores, 58 in the HCP and 36 in the ABCD).

### Running the code 
Create an output directory to store the images for the figures relating to HCP, and a directory for images for ABCD. Modify the results directory and the output directory for the images in part 3 of the code to suit your setup. The code then generates all figures in the output directory.
The code calls `HCP_variables_to_predict_real_names.txt` and `ABCD_variables_to_predict_real_names.txt`. These are the original behaviors in each dataset which has been expanded from their shortform. This code assumes that the order of these files and the order of the behaviors in the results file are the same.