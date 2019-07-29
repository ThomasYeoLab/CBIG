# Behaviour Partial Least Squares (PLS) for Neuroimaging

The main function is `myPLS_analysis.m`, which will run the following steps:

## PLS steps
- Data normalization
- Cross-covariance matrix
- Singular value decomposition
- Brain & behavior scores
- Brain & behavior loadings 
- Permutation testing for components' significance

## Requirements

Requires SPM8/12 for loading & reading volumes, and for visualisation (@slover functions).


## Credits

Code written by Dimitri Van De Ville, Daniela Zöller, Valeria Kebets and Thomas AW Bolton, with subfunctions borrowed from the PLS toolbox by Rotman Baycrest (https://www.rotman-baycrest.on.ca/index.php?section=84).

Please cite the following papers when using this code:

Zoller D, Schaer M, Scariati E, Padula MC, Eliez S, Van De Ville D (2017). Disentangling resting-state BOLD variability and PCC functional connectivity in 22q11.2 deletion syndrome. Neuroimage 149, pp. 85-97.

McIntosh AR, Lobaugh NJ (2004). Partial least squares analysis of neuroimaging data: applications and advances. Neuroimage 23(Suppl 1), pp. S250-263.	


If you have issues, please email Valeria Kebets (valkebets@gmail.com)
