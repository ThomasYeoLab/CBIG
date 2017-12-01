## Reference

Yeo BT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, Roffman JL, Smoller JW, Zollei L., Polimeni JR, Fischl B, Liu H, Buckner RL. [**The organization of the human cerebral cortex estimated by intrinsic functional connectivity**](http://www.ncbi.nlm.nih.gov/pubmed/21653723). Journal of Neurophysiology 106(3):1125-1165, 2011.

Krienen FM, Yeo BTT, Buckner RL. [**Reconfigurable state-dependent functional coupling modes cluster around a core functional architecture**](http://people.csail.mit.edu/ythomas/publications/2014TaskDependentCouplingModes-PTBS.pdf). Philosophical Transactions of the Royal Society B, 369:20130526, 2014.

Yeo BTT, Tandi J, Chee MWL. [**Functional connectivity during rested wakefulness predicts vulnerability to sleep deprivation**](http://people.csail.mit.edu/ythomas/publications/2015SleepDeprivation-NeuroImage.pdf). Neuroimage 111:147-158, 2015. 

----

## Updates
- Release v0.4.3 (09/10/2017): Release Yeo2011 brain parcellation.
- Release v0.4.5 (01/12/2017):

	**Yeo2011_fcMRI_clustering**
	
	1. Add Yeo2011 parcellations with split labels.
	
	**CBIG_fMRI_preprocessing**
	
	1. Change motion correction (mcflirt) interpolation method from default **trilinear** to **spline**.
	
	2. Add an optional preprocessing step to perform despiking by **AFNI 3dDespike**.
	
	3. Add a preprocessing step to generate ROIs2ROIs functional connectivity matrix for input subject. 
 
----

## Parcellation Release
The pacellations are released in `1000subjects_reference` folder. Specifically, the pacellations include:
- **1000subjects_clusters007_ref.mat, 1000subjects_clusters017_ref.mat**
7/17-network brain pacellation

- **Yeo_JNeurophysiol11_SplitLabels**
A connected component analysis was performed on the original 7/17-network brain parcellations. For more information, please check `Yeo_JNeurophysiol11_SplitLabels_README` under `Yeo_JNeurophysiol11_SplitLabels` folder

----

## Code Release
The code utilized in this study are released in `Yeo2011_fcMRI_clustering` folder. Specifically, the code include:
- **CBIG_general_cluster_fcMRI_surf2surf_profiles.csh**
This function is the main function that calls other scripts in sequence. It assumes that the preprocessed surface data are located in `sub_dir/subject/surf/`. Try `./CBIG_general_cluster_fcMRI_surf2surf_profiles.csh -help` for more information.

Note that this project uses generic functions from other folders, which may be updated over time. To download the version of the code that was last tested, you can either

- visit this link:
[https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.4.5-Yeo2011_fcMRI_clustering](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.4.5-Yeo2011_fcMRI_clustering)

or

- run the following command, if you have Git installed
 
```
git checkout -b v0.4.5-Yeo2011_fcMRI_clustering v0.4.5-Yeo2011_fcMRI_clustering
```

----

## Bugs and Questions

Please contact Thomas Yeo at yeoyeo02@gmail.com.

