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
For each parcel resolution there are two nifti files corresponding to different versions of MNI space. `FSLMNI152_1mm` and `FSLMNI152_2mm` contain parcellations which are aligned with 2mm and 1mm MNI templates coming FSL. These two FSL MNI templates can be found under `${CBIG_CODE_DIR}/data/templates/volume/`. 

The parcellations were computed in fsaverage6 space and sampled to the volumetric MNI space. Each parcel was matched to a corresponding network in the 7 and 17 network parcellation by Yeo et al. 

----

Example Usage
=============
### Freeview
1) Make sure freesurfer has been set up and configured as usual (http://surfer.nmr.mgh.harvard.edu/fswiki/SetupConfiguration).  
 
2) In terminal: 

   a) `cd` to unzipped folder containing this README file, Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz  
   
   b) enter the following command:

      ```
      freeview -v ${CBIG_CODE_DIR}/data/templates/volume/FSL5.0.8_MNI_templates/MNI152_T1_2mm_brain.nii.gz   ./freeview_lut/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz:colormap=lut:lut=Schaefer2018_400Parcels_7Networks_order.txt
      ```  
   
   c) set colormap from `Grayscale` to `Lookup Table`  


### fsleyes
1) Make sure FSL has been set up and configured as usual (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/ShellSetup)
 
2) In terminal: 

   a) `cd` to unzipped folder containing this README file, Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii.gz
  
   b) open `fsleyes`: `File -> Add from file -> Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii.gz`
  
   c) In the overlay information panel, set overlay type as `Label image`. By default, fsleyes will set it as `3D/4D volume`

   d) If your lookup tabels panel is not shown, go to `Settings -> Ortho View1 -> Lookup tabels`

   e) In the lookup tabels panel, select `Load LUT`, choose `fsleyes_lut/Schaefer2018_400Parcels_17Networks_order.lut`. Wait for a few seconds then you will see the colored parcellation.

Check for more details in the screenshot:
![Schaefer_fsleyes](https://user-images.githubusercontent.com/20438248/118851796-ba356000-b904-11eb-8f33-33d05cd7fa8a.png) 

----

Bugs and Questions
==================
Please contact Alexander Schaefer at alexschaefer83@gmail.com and Thomas Yeo at yeoyeo02@gmail.com.


