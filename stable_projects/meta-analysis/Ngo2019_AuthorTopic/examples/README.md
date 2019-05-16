# Example code for running author-topic model

---- 

## Reference
Gia H. Ngo, Simon B. Eickhoff, Minh Nguyen, Peter T. Fox,  R. Nathan Spreng, B. T. Thomas Yeo. [Beyond Consensus: Embracing Heterogeneity in Curated Neuroimaging Meta-Analysis](https://www.biorxiv.org/content/early/2017/06/19/149567). BioRxiv preprint

Gia H. Ngo, Simon B. Eickhoff, Peter T. Fox, B. T. Thomas Yeo. [Collapsed variational Bayesian inference of the author-topic model: application to large-scale coordinate-based meta-analysis](https://ieeexplore.ieee.org/document/7552332). PRNI2016.

---

## Overview
`CBIG_AuthorTopic_RunExample.m` is an example Matlab file to perform an end-to-end coordinate-based meta-analysis using the author-topic model, performing the following steps:

1. Preprocess MNI152 activation coordinates from a text file.
2. Estimate model parameters using Collapsed Variational Bayes inference.
3. Compute Bayesian Information Criterion (BIC) to determine the best number of components.
4. Generate surface-based visualization of the components.

`CBIG_AuthorTopic_RunExample.m` can be run out of the box. However, for actual meta-analytic experiment, larger range of number of components and higher number of re-initializations should be used to obtain stable estimates. Please refer to comments in `CBIG_AuthorTopic_RunExample.m` for the recommended numbers.

1. [SelfGeneratedThought](./SelfGeneratedThought) folder contains data and instruction to discover cognitive components of the functional domain of self-generated thought.
2. [CoactivationMappingIFJ](./CoactivationMappingIFJ) folder contains instruction to estimate the co-activation patterns of a seed region such as the inferior frontal junction (IFJ).

----

## Functional Sub-domains Discovery

`CBIG_AuthorTopic_RunExample.m` uses activation coordinates of 7 tasks engaging self-generated thought as input. To replicate the results we obtained, only two hyperparameters need to be changed:

```
  allKs = 1:4; % Line 50, set the number of components to be estimated to be between 1 and 5 inclusive
  seeds = 1:2; % Line 51, use 1000 reinitializations to obtain stable estimates
```

Please refer to [../SelfGeneratedThought/README.md](../SelfGeneratedThought) folder for details on the format of the input data.

----

## Co-activation Mapping

`CBIG_AuthorTopic_RunExample.m` can be applied to perform co-activation mapping without any major changes to the code (besides setting the number of components and reinitializations that you desire).

The only difference when applying the author-topic model to co-activation mapping as compared to functional sub-domains discovery is that: while experiments under a functional domain employing a limited number of behaviour tasks under the given domain (e.g. 179 experiments of self-generated thought employing 7 behavioral tasks), each of the experiments activating a seed region was considered to employ its own unique task category (e.g. 323 experiments activating the IFJ was considered as employing its own unique task category).

[../CoactivationMappingIFJ/README.md](../CoactivationMappingIFJ) describes the required input format for co-activation mapping meta-analysis, and instructions for replicating our experiment to estimate the co-activation patterns of the IFJ.
