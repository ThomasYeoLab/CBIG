# FSL-VBM -- Specific to Our PNAS Paper

This folder allows you to replicate the VBM results in our PNAS paper, provided that `../../extractBrains/replicatePNAS/` is used for skull stripping.

The wrapper first runs VBM from scratch on 810 ADNI-1 baseline scans and then runs VBM on 560 Month 24 (m24) follow-up scans using the study-specific template previously created with the 810 scans.

Another purpose of this folder is to provide examples of how to invoke `../lib/CBIG_runVBM.sh` and `../lib/CBIG_runVBM_givenTmp.sh`.

To learn the general procedure that can be easily applied to another dataset, please see folder `../lib/` instead.

