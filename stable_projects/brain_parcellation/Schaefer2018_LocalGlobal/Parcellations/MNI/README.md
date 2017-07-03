Background
==========
Resting state fMRI data from 1489 subjects were registered using surface-based alignment. A gradient weighted markov random field approach was employed to identify cortical areas.

Information about Downloads
===========================
For each parcel resolution there are three nifti files corresponding to different versions of MNI space. "FSLMNI152_1mm" and "FSLMNI152_2mm" contain parcellations which are aligned with 2mm and 1mm MNI templates coming FSL. These FSL MNI templates are also inside this folder with the names MNI152_T1_1mm_brain.nii.gz and MNI152_T1_2mm_brain.nii.gz. The parcellations were computed in "fsaverage6" space and sampled to the volumetric MNI space.  Each parcel was matched to a corresponding network in the 7 and 17 network parcellation by Yeo et al. 

Example Usage
=============

### Freeview
 1) Make sure freesurfer has been set up and configured as usual (http://surfer.nmr.mgh.harvard.edu/fswiki/SetupConfiguration).  
 2) In terminal  
   a) "cd" to unzipped folder containing this README file, Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz  
   b) ```freeview -v MNI152_T1_2mm_brain.nii.gz   Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz:colormap=lut:lut=Schaefer2018_400Parcels_7Networks_order.txt```  
 c) set colormap from 'Grayscale' to 'Lookup Table'  

### fslview
 1) Make sure FSL has been set up and configured as usual (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/ShellSetup)
 2) In terminal,  
  a) "cd" to unzipped folder containing this README file, Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz  
  b) ```fslview Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz -l Schaefer2018_400Parcels_7Networks_order.lut```  
  c) set min to 1 and max to 401  

References
==========
* Schaefer A, Kong R, Gordon EM, Zuo XN, Holmes AJ, Eickhoff SB, Yeo BT, (accepted), Local-Global Parcellation of the Human Cerebral Cortex From Intrinsic Functional Connectivity MRI, Cerebral Cortex
* Yeo BT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, Roffman JL, Smoller JW, Zollei L., Polimeni JR, Fischl B, Liu H, Buckner RL (2011) The organization of the human cerebral cortex estimated by functional connectivity. J. Neurophysiol.
