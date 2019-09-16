References
==========
+ Schaefer A, Kong R, Gordon EM, Laumann TO, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BTT. [**Local-Global parcellation of the human cerebral cortex from intrinsic functional connectivity MRI**](http://people.csail.mit.edu/ythomas/publications/2018LocalGlobal-CerebCor.pdf), *Cerebral Cortex*, 29:3095-3114, 2018

----

Parcellation Code
===========================
The code in this folder was used to project the parcellations from ```fsaverage``` space to 
```MNI152``` and ```fslr32k``` space. The parcellations in ```fsaverage``` space are contained in the ```input``` subfolder. The code in the current form does only work with the data in the ```input``` folder.

----

Usage
=============
To assign parcel labelnames to Schaefer parcellations, Our matching algorithm does the matching in two steps:

1. The Schaefer parcels are first matched to one of Yeo2011 7/17 networks. Each parcel from Schaefer parcellations is then assigned to a network, e.g., `17Networks_LH_SomMotB`. 
2. Each Schaefer's parcel is then matched to one of the Yeo2011 split components and is assigned a component name, e.g., `17Networks_LH_SomMotB_Cent`.

Function `CBIG_gwMRF_regenerate_Schaefer2018_parcellations.m` reads in the Schaefer parcellations from current repo and performs the matching procedure to generate a new version of Schaefer2018 parcellations. The parcel labelnames and the ordering of Schaefer parcellations might change in the new version. 

Run the following command to regenerate Schaefer2018 parcellations under given folder:

```
CBIG_gwMRF_regenerate_Schaefer2018_parcellations('example_folder')
```

To copy common label and surface files to the folder, a script can be used: 
```
sh CBIG_gwMRF_copy_fs_average.sh
```

----

Mapping parcellation ordering between versions
=============
Since the parcellation indices might change between different releases, we provide a function that compares two versions and generates a vector of indices to map the old parcellation ordering to the new parcellation ordering. Please refer to `CBIG_gwMRF_index_trans_btwn2versions.m`

----

Bugs and Questions
==================
Please contact Alexander Schaefer at alexschaefer83@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

