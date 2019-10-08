# Behaviour Partial Least Squares (PLS) for Neuroimaging

The main function is `myPLS_analysis.m`, which will run the following steps:

## PLS steps
- Data normalization
- Cross-covariance matrix
- Singular value decomposition
- RSFC & behavior composite scores
- RSFC & behavior loadings 

## Requirements
Requires SPM8/12 for loading & reading volumes, and for visualisation (@slover functions).


## Credits
These functions were modified from the original code written by Dimitri Van De Ville, Daniela Zöller, and Valeria Kebets.
The code is maintained by Daniela Zöller and Valeria Kebets, and is available at https://github.com/danizoeller/myPLS). In addition to myPLS functions, it includes subfunctions borrowed from the PLS toolbox by Rotman Baycrest (https://www.rotman-baycrest.on.ca/index.php?section=84).

If you have issues, please email Valeria Kebets (valkebets@gmail.com).
