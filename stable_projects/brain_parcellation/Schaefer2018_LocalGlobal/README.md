References
==========
+ Schaefer A, Kong R, Gordon EM, Laumann TO, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BTT. [**Local-Global parcellation of the human cerebral cortex from intrinsic functional connectivity MRI**](http://people.csail.mit.edu/ythomas/publications/2018LocalGlobal-CerebCor.pdf), *Cerebral Cortex*, 29:3095-3114, 2018

----

Background
==========
Resting state fMRI data from 1489 subjects were registered using surface-based alignment. A gradient weighted markov random field approach was employed to identify cortical parcels ranging from 100 to 1000 parcels.
More details can be found in Schaefer et al. 2018.

---

Parcellations Release
=====================
The parcellations are available at multiple resolution (100 parcels to 1000 parcels), and can be found under the ```Parcellations``` folder. To use the parcellations without the trouble of downloading our entire repository, you can just click on this link: [download Schaefer2018_Parcellations](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations)

Specifically, there are three subfolders corresponding to three different spaces ```Freesurfer5.3```, ```MNI``` and ```HCP```. 

The parcellations were first computed in Freesurfer ```fsaverage6``` space and sampled to ```fsaverage5``` and ```fsaverage``` space. The parcellations were also projected to HCP ```fslr32k``` and FSL ```MNI``` space. Each parcel is matched to a corresponding network in the 7 and 17 network parcellation by Yeo et al. 2011.([The organization of the human cerebral cortex estimated by intrinsic functional connectivity](http://www.ncbi.nlm.nih.gov/pubmed/21653723)).

Here we provide a visualization of the 400 parcel parcellation in ```fslr32k``` space, parcels were colored to match Yeo 7/17 network parcellation:

<img src="readme_figures/Schaefer2018_400parcel_parcellation_match_Yeo_7_network_fslr32k.png" height="300" />

<img src="readme_figures/Schaefer2018_400parcel_parcellation_match_Yeo_17_network_fslr32k.png" height="300" />


We also provide RAS centroid coordinates of the parcellations in MNI 1mm and 2mm space. If you are interested, please check the csv files in: `Parcellations/MNI/Centroid_coordinates`

---

Code Release
============
The code utilized in this study is under `Code` folder. Specifically, the `Code` folder includes:

* **CBIG_gwMRF_build_data_and_perform_clustering.m** -- The wrapper function that generates input data format from fMRI data and perform clustering. 

* **lib** folder -- This folder contains all other functions that will be called by the wrapper function.

* **README.md** -- You can check this file to find out more about our clustering code.

### Examples
We provide example code for you run as well as example results for you to compare under `examples` folder. Please refer to `examples/README.md` for more information.

### Download stand-alone repository
Since the whole Github repository is too big, we provide a stand-alone version of only this project and its dependencies. To download this stand-alone repository, visit this link: [https://github.com/ThomasYeoLab/Standalone_Schaefer2018_LocalGlobal](https://github.com/ThomasYeoLab/Standalone_Schaefer2018_LocalGlobal)

### Download whole repository
Except for this project, if you want to use the code for other stable projects from out lab as well, you need to download the whole repository.

To download the version of the code that is last tested, you can either

- visit this link: [https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.8.1-Schaefer2018_LocalGlobal](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.8.1-Schaefer2018_LocalGlobal)

or

- run the following command, if you have Git installed

```
git checkout -b Schaefer2018_LocalGlobal v0.8.1-Schaefer2018_LocalGlobal
```

---

Updates
=======

- Release v0.3.0 (21/07/2017): release the parcellations in MNI, fsLR and fsaverage space

- Release v0.3.1 (07/09/2017): update utilities code

- Release v0.3.2 (14/09/2017): fix a labeling bug in cifti files

- Release v0.4.8 (30/01/2018): release clustering code

- Release v0.4.12 (09/04/2018): 

    1. Move example subjects from `$CBIG_CODE_DIR/data/example_data/${subj_ID}` to `$CBIG_CODE_DIR/data/example_data/Corr_HNU/${subj_ID}`. 
    
    2. The example subjects are re-processed by a newer version of our preprocessing pipeline ([v0.4.9](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.4.9-CBIG_fMRI_Preprocessing)). Hence results in `./Code/examples` are updated.
    
    3. Update some README.md files.

- Release v0.6.3 (17/07/2018):

    1. Release Schaefer2018 300 parcels and 500 parcels parcellation.

    2. Add unit test scripts to `unit_tests` folder.

    3. Update README.md, add links for people to directly download the Parcellations and this project's standalone repository: `Standalone_Schaefer2018_LocalGlobal`.

- Release v0.6.5 (17/08/2018):

    1. Move `examples` folder from `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code/examples` to `$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/examples`.

    2. Update publication details in README files.

    3. Fix a bug related to creating `Standalone_Schaefer2018_LocalGlobal`.

- Release v0.8.1 (01/02/2019):

    1. Release RAS centroid coordinates of Schaefer2018 parcellations in MNI 1mm and 2mm space. The coordinate csv files can be found under `Parcellations/MNI/Centroid_coordinates`.

---

Bugs and Questions
==================
Please contact Alexander Schaefer at alexschaefer83@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

