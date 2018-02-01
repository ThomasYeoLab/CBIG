References
==========
+ Schaefer A, Kong R, Gordon EM, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BT, (accepted), Local-Global Parcellation of the Human Cerebral Cortex From Intrinsic Functional Connectivity MRI, Cerebral Cortex

----

Background
==========
Resting state fMRI data from 1489 subjects were registered using surface-based alignment. A gradient weighted markov random field approach was employed to identify cortical parcels.
More details can be found in Schaefer et al. 2018.

---

Updates
=======
- Release v0.4.8 (30/01/2018): release clustering code

- Release v0.3.2 (14/09/2017): fix a labeling bug in cifti files

- Release v0.3.1 (07/09/2017): update utilities code

- Release v0.3.0 (21/07/2017): release the parcellations in MNI, fsLR and fsaverage space

---

Parcellations Release
=====================
The parcellations can be found in ```Parcellations``` folder. There are three subfolders corresponding to three different spaces ```Freesurfer5.3```, ```MNI``` and ```HCP```. 
The parcellations were computed in Freesurfer ```fsaverage6``` space and projected to HCP ```fslr32k``` and FSL ```MNI``` space. Each parcel is matched to a corresponding network in the 7 and 17 network parcellation by Yeo et al. 2011.  

---

Code Release
============
The code utilized in this study is under `Code` folder. Specifically, the `Code` folder includes:

* **CBIG_gwMRF_build_data_and_perform_clustering.m** -- The wrapper function that generates input data format from fMRI data and perform clustering. 

* **lib** folder -- This folder contains all other functions that will be called by the wrapper function.

* **examples** folder -- In this folder we provide example code for you to run and example results for you to compare. Please refer to `/Code/examples/README.md` for more information.

* **README.md** -- You can check this file to find out more about our clustering code.


To download the version of the code that is last tested, you can either

- visit this link: [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.4.8-Schaefer2018_LocalGlobal](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.4.8-Schaefer2018_LocalGlobal)

or

- run the following command, if you have Git installed

```
git checkout -b Schaefer2018_LocalGlobal v0.4.8-Schaefer2018_LocalGlobal
```

---

Bugs and Questions
==================
Please contact Alexander Schaefer at alexschaefer83@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

