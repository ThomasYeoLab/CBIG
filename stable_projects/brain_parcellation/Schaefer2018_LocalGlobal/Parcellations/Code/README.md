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
Code can be exectuted calling:
```
CBIG_gwMRF_run_Schaefer2018_create_all_public_parcellations('example_folder/')
```

To copy common label and surface files to the folder, a script can be used: 
```
sh CBIG_gwMRF_copy_fs_average.sh
```

----

Bugs and Questions
==================
Please contact Alexander Schaefer at alexschaefer83@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

