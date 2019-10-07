This folder compares RSFC between all ASD and control participants, as well as RSFC between ASD and control participants in each factor subgroup in ABIDE-I. Statistical significance is tested using network-based statistic (NBS; Zalesky et al. 2010).

Figures S5 & S6 are relevant to this folder.

----
## What Does Each File Do?
1. `CBIG_ASDf_genPermSetForNBS.m` generates within-site permutation set using PALM package.
2. `CBIG_ASDf_NBS_permSetInput.m` performs NBS taking in permutation set as input to allow for within-site permutation.
3. `CBIG_ASDf_FCDiffAllSub_NBS.m` compares RSFC difference between all ASD & control participants. Statistical significance is tested using NBS.
4. `CBIG_ASDf_NBSall_job.sh` calls `CBIG_ASDf_FCDiffAllSub_NBS.m` and submits the job to our cluster.
5. `CBIG_ASDf_FCDiffInSubgrp_NBS.m` compares RSFC difference between ASD & control participants in factor subgroup. Statistical significance is tested using NBS.
6. `CBIG_ASDf_NBSsubgrp_job.sh` calls `CBIG_ASDf_FCDiffInSubgrp_NBS.m` and submits the job to our cluster.
7. `CBIG_ASDf_plotFCDiff_NBSThresholded.m` plots RSFC difference between ASD & control participants after performing NBS.
8. `CBIG_ASDf_FCDiffInFactorGroups_wrapper.m` is the wrapper function to compare RSFC between all ASD and control participants, as well as RSFC between ASD and control participants in each factor subgroup in ABIDE-I.
