This folder generates Table S1.

This folder investigates, with volume estimations by FreeSurfer, whether an atrophy factor is indeed associated with severe atrophy in the corresponding brain regions.

---

## What Does Each File Do?

1. `CBIG_wrapper.m` calls all other functions here.
2. `CBIG_getVol.m` fetches brain structure volumes from FreeSurfer `stats` files.
3. `CBIG_sortByAvgProbWinningFactor.m` sorts structure names for displaying in the table.
4. `CBIG_assignFactorsToStructures.m` assigns each brain structure to an atrophy factor.
5. `CBIG_fitGLM_hypoTest.m` fits a GLM and then tests whether atrophy factors really reflect atrophy in the respective brain regions.
