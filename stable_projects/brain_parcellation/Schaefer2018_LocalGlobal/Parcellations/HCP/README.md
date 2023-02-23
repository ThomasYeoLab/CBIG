References
==========
+ Schaefer A, Kong R, Gordon EM, Laumann TO, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BTT. [**Local-Global parcellation of the human cerebral cortex from intrinsic functional connectivity MRI**](http://people.csail.mit.edu/ythomas/publications/2018LocalGlobal-CerebCor.pdf), *Cerebral Cortex*, 29:3095-3114, 2018
+ Yeo BTT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, Roffman JL, Smoller JW, Zöllei L, Polimeni JR, Fischl B, Liu H, Buckner RL. [**The organization of the human cerebral cortex revealed by intrinsic functional connectivity**](http://people.csail.mit.edu/ythomas/publications/2011CorticalOrganization-JNeurophysiol.pdf), *J Neurophysiology*, 106(3):1125-1165, 2011


----

Background
==========
Resting state fMRI data from 1489 subjects were registered using surface-based alignment. A gradient weighted markov random field approach was employed to identify cortical areas.

----

Information about Downloads
===========================
For each parcel resolution there are two `dlabel.nii` files corresponding to assignments to Yeo et al 7 and 17 network parcellation. The parcellations were computed in "fsaverage6" space and sampled to the fsLR32k space.  Each parcel was matched to a corresponding network in the 7 and 17 network parcellation by Yeo et al. 

----

Example Usage
=============
### Workbench Viewer
1) Make sure connectome workbench viewer has been set up and configured (https://wiki.humanconnectome.org/display/WBPublic/Connectome+Workbench+Help+Guide).

2) In terminal:

   a) `wb_view`
   
   b) load fslr32k surface files for both hemispheres. For example, you can use these files under this directory `CBIG/data/templates/surface/fs_LR_32k/` in our [public repo](https://github.com/ThomasYeoLab/CBIG):
   - [fsaverage.L.very_inflated.32k_fs_LR.surf.gii](https://github.com/ThomasYeoLab/CBIG/blob/master/data/templates/surface/fs_LR_32k/fsaverage.L.very_inflated.32k_fs_LR.surf.gii)
   - [fsaverage.R.very_inflated.32k_fs_LR.surf.gii](https://github.com/ThomasYeoLab/CBIG/blob/master/data/templates/surface/fs_LR_32k/fsaverage.R.very_inflated.32k_fs_LR.surf.gii)
   
   c) select file-> open file  
   
   d) click to the folder containing this Readme and from there click to subfolders `fslr32k` and `cifti` containing the `dlabel.nii` files
   
   e) change `Files of type` to `Connectivity Dense Label Files (*dlabel.nii)`
   
   f) select the Schaefer_2018 parcellation file you want to use  
   
   g) in the overlay toolbox select the corresponding Schaefer_2018 parcellation file


### Matlab
1) Make sure Matlab has been installed together with the FieldTrip Cifti Reader (https://github.com/Washington-University/cifti-matlab)  

2) Run the following command:
   ```
   x = ft_read_cifti('Schaefer2018_400Parcels_17Networks_order.dlabel.nii','mapname','array');
   ```

----

Bugs and Questions
==================
Please contact Alexander Schaefer at alexschaefer83@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

