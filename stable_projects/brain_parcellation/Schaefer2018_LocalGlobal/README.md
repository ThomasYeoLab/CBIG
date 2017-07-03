Background
==========
Resting state fMRI data from 1489 subjects were registered using surface-based alignment. 
A gradient weighted markov random field approach was employed to identify cortical parcels.
More details can be found in Schaefer et al. 2018.

Parcellations Downloads
===========================
The parcellations can be found in ```Parcellations```. There are three subfolders corresponding to three different 
spaces ```Freesurfer5.3```, ```MNI``` and ```HCP```. The parcellations were computed in Freesurfer ```fsaverage6``` space and projected to 
HCP ```fslr32k``` and FSL ```MNI``` space. Each parcel was matched to a corresponding network in the 7 and 17 network parcellation by Yeo et al. 2011.  

Release
=======
- visit this link:
[https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.3.0-Schaefer2018_LocalGlobal](https://github.com/ThomasYeoLab/CBIG/releases/tag/v0.3.0-Schaefer2018_LocalGlobal)

or

- run the following command, if you have Git installed

```
git checkout -b Schaefer2018_LocalGlobal v0.3.0-Schaefer2018_LocalGlobal
```

Bugs and Questions
==================
Please contact Alexander Schaefer at alexschaefer83@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

References
==========
+ Schaefer A, Kong R, Gordon EM, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BT, (accepted), Local-Global Parcellation of the Human Cerebral Cortex From Intrinsic Functional Connectivity MRI, Cerebral Cortex
