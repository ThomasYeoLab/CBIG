# Disclaimer
This is the pre-alpha version of the code and data to reproduce the experiments reported in our preprint. We are working hard on a stable and user-friendly release that can be easy to run on new datasets.

# Reference
Gia H. Ngo,  Simon B. Eickhoff, Peter T. Fox,  R. Nathan Spreng, B. T. Thomas Yeo. [Beyond Consensus: Embracing Heterogeneity in Neuroimaging Meta-Analysis](https://app.bitly.com/Bg84gvJLXEM/bitlinks/2ruE7v8). BioRxiv preprint

## Data Release 
The two sets of activation coordinates used in our experiments were included under the folder `self-generated_thought/input_data` and `IFJ/input_data`.


## Code Release
1. `utilities` contains all the functions necessary to replicate our exeperiments from scracth: estimating the cognitive components/ co-activation patterns, compute Bayesian Information Criterion and visualize the results.
2. `self-generated_thought` contains the wrapper scripts necessary to reproduce the cognitive components of self-generated thought.
3. `IFJ` contains the wrapper scritps necessary to reproduce the co-activation patterns that the IFJ expresses.

Note that this project uses generic functions from other folders, which may be updated over time. To download the exact code used for our experiments, you can either

- visit this link:
[https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.2.0-Ngo2017_AuthorTopic](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.2.0-Ngo2017_AuthorTopic)

or

- run the following command, if you have Git installed
 
```
git checkout -b Ngo2017_AuthorTopic v0.2.0-Ngo2017_AuthorTopic
```

----

## Bugs and Questions

Please contact Gia at ngohoanggia@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.
