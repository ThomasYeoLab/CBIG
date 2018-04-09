This folder contains subjects from the preprocessed CoRR_HNU dataset. We use these subjects as example data so that people can run our example code and make sure the outputs are correct.

----

CoRR_HNU dataset
================
The CoRR_HNU dataset is an open source science resource created by CoRR (Consortium for Reliability and Reproducibility). We get permission from Professor Xi-Nian Zuo to re-release a few preprocessed CoRR_HNU subjects here as example data.

You may refer to this paper for more information about the CoRR_HNU dataset: [An open science resource for establishing reliability and reproducibility in functional connectomics](https://www.nature.com/articles/sdata201449.pdf).

----

CBIG preprocessing pipeline
===========================
Example data from the CoRR_HNU dataset have been preprocessed with our CBIG pipeline (version number v0.4.9).

The order of preprocessing steps are listed below:

* CBIG_preproc_skip -skip 4

* CBIG_preproc_fslslicetimer -slice_order ${CBIG_CODE_DIR}/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/example_slice_order.txt

* CBIG_preproc_fslmcflirt_outliers -FD_th 0.2 -DV_th 50 -discard-run 50 -rm-seg 5 -spline_final

* CBIG_preproc_bbregister -intrasub_best

* CBIG_preproc_regress -whole_brain -wm -csf -motion12_itamar -detrend_method detrend -per_run -censor -polynomial_fit 1

* CBIG_preproc_censor -max_mem NONE

* CBIG_preproc_bandpass -low_f 0.009 -high_f 0.08 -detrend 

* CBIG_preproc_QC_greyplot -FD_th 0.2 -DV_th 50

* CBIG_preproc_native2fsaverage -proj fsaverage6 -sm 6 -down fsaverage5


