Information about Downloads
===========================
The parcellations were computed in "fsaverage6" space and sampled to the fsLR32k space.


Example Usage
=============
### Workbench Viewer
1) Make sure connectome workbench viewer has been set up and configured (https://wiki.humanconnectome.org/display/WBPublic/Connectome+Workbench+Help+Guide).

2) In terminal:

   a) `wb_view`  
   
   b) load fslr32k surface files for both hemispheres. For example, you can use these files under this directory `CBIG/data/templates/surface/fs_LR_32k/` in our [public repo](https://github.com/ThomasYeoLab/CBIG):
   - [fsaverage.L.very_inflated.32k_fs_LR.surf.gii](https://github.com/ThomasYeoLab/CBIG/blob/master/data/templates/surface/fs_LR_32k/fsaverage.L.very_inflated.32k_fs_LR.surf.gii)
   - [fsaverage.R.very_inflated.32k_fs_LR.surf.gii](https://github.com/ThomasYeoLab/CBIG/blob/master/data/templates/surface/fs_LR_32k/fsaverage.R.very_inflated.32k_fs_LR.surf.gii)
   
   c) select file -> open file  
   
   d) click to the folder containing this Readme and from there click to subfolders `fslr32k` containing the `dlabel.nii` files
   
   e) change `Files of type` to `Connectivity Dense Label Files (*dlabel.nii)`
   f) select the parcellation file you want to use  
   
   g) in the overlay toolbox, select the loaded parcellation file


### Matlab
1) Make sure Matlab has been installed together with the FieldTrip Cifti Reader (https://github.com/Washington-University/cifti-matlab)  

2) Run the following command:
   ```
   x = ft_read_cifti('100Parcels_Kong2022_17Networks.dlabel.nii','mapname','array');
   ```

Bugs and Questions
==================
Please contact Xiaoxuan Yan at xiaoxuan.427@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.

