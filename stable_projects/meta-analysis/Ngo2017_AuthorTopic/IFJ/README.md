This folder contains wrapper functions necessary to reproduce the results on co-activation patterns that the inferior frontal junction (IFJ) expresses (Fig. 5)

# Initialization
From Matlab, run:
```
initCVB; % add the necessary paths
```

# Inference
To peform inference with the number of co-activation patterns K ranging from 1 to 4, with 100 random initialization for each K, run in Matlab the following command:
```
CBIG_RunInferenceWrapper(1:4, 1:100);
```
The 100 estimates from the 100 initializations for each value of K would be stored under `outputs/K1` to `outputs/K4`

# Get the best estimate
To get the estimate with the highest variational bound from the 100 estimates for each value of K ranging from 1 to 4, run in Matlab the following command from this directory:
```
CBIG_ComputeBestSolutionWrapper(1:4, 1:100);
```
The best solution for each value of K would be stored under `outputs/best_solution`

# Compute Bayesian Information Criterion (BIC)
To compute the BIC value corresponding to each value of K ranging from 1 to 4, run the followign command in Matlab from this directory:
```
CBIG_ComputeBICWrapper(1:4)
```
This would produce the BIC plot in Fig. 4a.

# Visualize the co-activation patterns:
To project the co-activation patterns to the brain surface for visualization, run the following commands in Matlab from this directory:
```
K = 3;
CBIG_VisualizeComponentsWrapper(K);
```
This would produce the surface visualization of the 3 co-activation patterns in Fig. 4b. The images are available under `outputs/best_solution/alpha100_eta0.01/figures`. The brain volumes corresponding to the 3 co-activation patterns are available udner `outputs/BIC/ComponentsSmoothness/`

# Get task loading of each co-activation pattern
The probability of each experiment recruiting a co-activation pattern is available in the `theta` matrix of the struct `params` stored in `outputs/best_solution/K2/alpha100_eta0.01/BestSolution_K003.mat`.
The 323 rows of the `theta` matrix correspond to 323 experiments that activates the inferior frontal junction from The BrainMap database.
See folder `input_data` for the 323 experiments.

