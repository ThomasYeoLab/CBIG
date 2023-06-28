References
==========
+ Schaefer A, Kong R, Gordon EM, Laumann TO, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BTT. [**Local-Global parcellation of the human cerebral cortex from intrinsic functional connectivity MRI**](http://people.csail.mit.edu/ythomas/publications/2018LocalGlobal-CerebCor.pdf), *Cerebral Cortex*, 29:3095-3114, 2018

----

Install
=======
**Set up your enviroment**

The configuration scripts `CBIG_gwMRF_tested_config.sh` and `CBIG_gwMRF_tested_startup.m` can be found under `../config` folder. 

Please follow the instructions in `<your-cbig-repo-direcotry>/setup/README.md` to have your local environment compatible with CBIG repository. 

----

Code
====
Schaefer2018 clustering code consists of a wrapper function and several lib functions it calls. Here we will briefly explain some core funtions:
1) Wrapper function:
```
CBIG_gwMRF_build_data_and_perform_clustering.m
```
This wrapper function generates input data format from fMRI data and perform clustering.


2) Functions under `lib` folder:
```
CBIG_gwMRF_build_time_matrix.m
```

This function concatenates your own input timeseries data.

```
CBIG_gwMRF_build_cov_matrix.m
```

This function premutiplies the timeseries data. This will reduce the memory needed if the timeseries are long. Check [Schaefer2018 supplementary](https://academic.oup.com/cercor/advance-article/doi/10.1093/cercor/bhx179/3978804?searchresult=1) for more details.

```
CBIG_gwMRF_set_prams.m
```

This function sets all the parameters used for clustering. 

```
CBIG_gwMRF_graph_cut_clustering_split_newkappa_prod.m
```

This function performs the actual clustering for the premulitiplied input data.

----  

Parameters
==========
1) Input arguments for the wrapper function `CBIG_gwMRF_build_data_and_perform_clustering.m`

  * `input_fullpaths`=  a file containing full paths to all subjects' surf data, each line represents a subject with different runs;
  * `output_path`=      path to where output files will be written;
  * `start_idx`=        for selecting subsets in the subject list, this will be the start index;
  * `end_idx`=          for selecting subsets in the subject list, this will be the end index;
  * `num_left_cluster`= number of cluster on the left hemisphere;
  * `num_righ_cluster`= number of cluster for the right hemisphere;
  * `smoothcost`=       weight for the smoothness in the MRF, the higher this value the more the gradient prior will be weighted in the optimization, see the paper for more details;
  * `num_iterations`=   number of iterations for each random initializations, it might stop earlier if the energy change is too small;
  * `num_runs`=         number of random initializations;
  * `start_gamma`=      this is the initial gamma, a roundness prior, the higher this value is the more round the parcells will be;
  * `exponential`=      this controls the steepnes of the gradient prior function (exp(-k * gradient) - exp(-k)). The higher the value the less stringent the gradient prior map will be.

2) Parameters defined in `CBIG_gwMRF_set_prams.m`

You can set many parameters for the clustering process in `CBIG_gwMRF_set_prams.m`. All parameters and their default setting should be explained in this function. 

Here we will explain some of the most commonly changed parameters, the rest you should rather leave as default or read in `CBIG_gwMRF_set_prams.m`.

  * `dim`=              the dimensionality of the original timeseries data, needed for premultiplicated input data;
  * `lh_avg_file`=      input data for the left hemisphere;
  * `rh_avg_file`=      input data for the right hemisphere;
  * `pca`=              1 if you want to use premultiplied data; (the term pca is misleading, should be named premultiplied)
  * `first_gamma`=      the first value gamma will be increased to; (to enforce a certain number of parcels, gamma might need to be increased)
  * `reduce_gamma`=     1 if you want to iteratively decrease gamma;
  * `reduce_speed`=     how fast the gamma will be reduced, in this case gamma'=gamma * 1/5.

----

Output
======
The clustering results will be under `<output_path>/clustering`, each random seed will create two output mat files in this folder.

For example for seed 1:
  1) `*_seed_1_Energy.mat`, contains only the energy (optimization result) and can be quickly loaded.

  2) `*_seed_1.mat`, contains a `prams` and `results` variable.
     * prams: contains all the parameters we defined above and are default
     * results: contains lh_label and rh_label which is our result and some other information (i.e. likelihood, Energy and so on).

----

Examples
========
We provide a simple example under `../examples` folder. You can check how to run the example code and compare the results by reading `../examples/README.md`.

----

Graph Cut optimizer
===================
We use the graph cut optimizer by Delong et al. from http://vision.csd.uwo.ca/code/ named gco-v3.0. 
More information can be found in the following paper "Fast Approximate Energy Minimization with Label Costs" International Journal of Computer Vision, vol. 96, no. 1, pp. 1–27, January 2012".

----

Bugs and Questions
==================
Please contact Alexander Schaefer at alexschaefer83@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.



