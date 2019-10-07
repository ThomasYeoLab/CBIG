This folder performs k-means clustering on z-normalized RSFC of ASD participants in ABIDE-II+GENDAAR, and compares behavioral associations in k-means clusters with that in latent factors.

Table S9 & Figure S4 are relevant to this folder.

----
## What Does Each File Do?
1. `CBIG_ASDf_clusterKmeans_wrapper.m` performs k-means clustering on the z-normalized RSFC of ASD participants in ABIDE-II+GENDAAR.
2. `CBIG_ASDf_plotKmeans.m` plots the cluster centroids on a 419x419 matrix.
2. `CBIG_ASDf_logReg_clusterBehavAssoc.m` and `CBIG_ASDf_clusterBehavAssoc_wrapper.m` perform logistic regression to find behavioral associations between k-means clusters and behavioral scores (Figure S4).
