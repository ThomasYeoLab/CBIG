Information about Downloads
===========================
For each parcel resolution there are two nifti files corresponding to different versions of MNI space. `FSLMNI152_1mm` and `FSLMNI152_2mm` contain parcellations which are aligned with 2mm and 1mm MNI templates.

Example Usage
=============
### Freeview
1) Make sure freesurfer has been set up and configured as usual (http://surfer.nmr.mgh.harvard.edu/fswiki/SetupConfiguration).  
 
2) In terminal: 

   a) Say now we want to visualize `100Parcels_Kong2022_17Networks_FSLMNI152_1mm.nii.gz`. First `cd` to the parent folder of this file.
   
   b) Enter the following command:

      ```
      freeview -v ${CBIG_CODE_DIR}/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_2mm_brain.nii.gz   100Parcels_Kong2022_17Networks_FSLMNI152_1mm.nii.gz:colormap=lut:lut=./freeview_lut/100Parcels_Kong2022_17Networks_LUT.txt
      ```  

### fsleyes
1) Make sure FSL has been set up and configured as usual (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/ShellSetup)
 
2) In terminal: 

   a) Say now we want to visualize `100Parcels_Kong2022_17Networks_FSLMNI152_1mm.nii.gz`. First `cd` to the parent folder of this file.
   
   b) Open `fsleyes`: `File -> Add from file -> 100Parcels_Kong2022_17Networks_FSLMNI152_1mm.nii.gz`
  
   c) In the overlay information panel, by default, the overlay is set to `3D/4D volume`. Toggle it to `Label image` instead.

   d) If your lookup tabels panel is not shown, go to `Settings -> Ortho View1 -> Lookup tabels`

   e) In the lookup tabels panel, select `Load LUT`, choose `fsleyes_lut/Schaefer2018_400Parcels_17Networks_order.lut`.

Check for an example in the screenshot:
![fsleyes](https://user-images.githubusercontent.com/20438248/118851796-ba356000-b904-11eb-8f33-33d05cd7fa8a.png) 


Bugs and Questions
==================
Please contact Xiaoxuan Yan at xiaoxuan.427@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.


