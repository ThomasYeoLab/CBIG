# Our LDA Files

These files are our model files that are necessary for inferring new subjects. 

---
## Model Files in `model_K*/`

These folders contain our models for K = 2, 3 and 4.

In each folder, `final.*` and `likelihood.dat` are the LDA outputs: `final.beta`, after exponentiated, gives the probability of a word (voxel) given a topic (factor); `final.gamma`, after normalization, gives the probability of a topic (factor) given a document (subject); `final.other` contains number of topics you asked for, number of words (voxels), and a model parameter after optimization; `likelihood.dat` shows how the likelihood evolves as inference proceeds.

For visualization purposes, we processed `final.beta` into brain volumes `factor*.nii.gz`, which you can load into FreeView or FSLView along with the MNI152 template `MNI152_T1_2mm_brain.nii.gz`. After performing inference, you will need to normalize `final.gamma` to obtain factor compositions that sum to 1 (such as `[0.7, 0.2, 0.1]`). To see how to process `final.beta` and `final.gamma`, see `../../../step2_LDA/lib/CBIG_visualizeFactors.m`.

----
## Parameters of the Reference Group `refParams.mat`

This `.mat` file contains regression and normalization parameters of the reference group.

Specifically, during the brain-to-document conversion, nuisance variables will be regressed out using a general linear model (GLM) whose parameters were estimated with the reference group (in our case, the 228 ADNI-1 CN participants). Additionally, the mean and standard deviations of the reference group are needed for the z-normalization. By loading this `.mat` file into MATLAB, you will have these necessary parameters to perform brain-to-document conversion.

The relevant code that uses these parameters is `../../../step2_LDA/lib/CBIG_brain2doc.m`.
