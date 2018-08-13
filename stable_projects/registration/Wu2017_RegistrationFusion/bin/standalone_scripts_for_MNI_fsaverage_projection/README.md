## Reference

Jianxiao Wu, Gia H. Ngo, Alexander Schaefer, Douglas Greve, Jingwei Li, Tong He, Bruce Fischl, Simon B. Eickhoff, B.T. Thomas Yeo. [**Accurate Nonlinear Mapping between MNI152/Colin27 Volumetric and FreeSurfer Surface Coordinate Systems**](http://people.csail.mit.edu/ythomas/publications/2018VolSurfMapping-HBM.pdf), *Human Brain Mapping*, 2018.

---

## Code Release

This folder contains standalone codes for using the Registration Fusion mappings for MNI-to-fsaverage projections. The mappings are produced in the same way as those used in our paper, but using all 1490 GSP subjects.

To download the codes for replication of our paper, or for creating customised mappings, you would need to install our Github repository (Refer to https://github.com/ThomasYeoLab/CBIG). You can find the project under stable_projects/registration/Wu2017_RegistrationFusion.

---

## Before you start

Make sure that FreeSurfer and Matlab are already installed. In addition, please make sure that you have set up FreeSurfer properly (see https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall#Setup.26Configuration for instructions).

 To test if the scripts work, you can try the example (see Example section).

---

## Example

Run `CBIG_RF_MNI_example.sh` for an example using `CBIG_RF_projectMNI2fsaverage.sh`. This projects a probabilistic map of central sulcus from MNI152 to fsaverage.

(Note that you may need to manually provide your Matlab path to the script using -m option. Use -h option for more information)

- this will generate two output files called `lh.MNI_probMap_ants.central_sulc.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii.gz` and `rh.MNI_probMap_ants.central_sulc.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii.gz`
- the lh output is compared to `lh.MNI_probMap_ants.central_sulc.demo.nii.gz` in Freeview. Check out the 2 files in `Overlay`. These 2 files should be the same, like shown below.

<p align="center">
<img src="stand_alone_MNI_example.png" height="300" />
</p>

The commands to run the same example are also included in the help message of the `CBIG_RF_projectMNI2fsaverage.sh` and `CBIG_RF_projectMNI2fsaverage.m`. 

--- 

## General Usage

To project data (A.nii.gz) in MNI152 space to fsaverage and save the outputs directly to folder (B_dir), use the Bash script with the command:

```bash
./CBIG_RF_projectMNI2fsaverage.sh -s full/path/to/A.nii.gz -o full/path/to/B_dir
```

(Note that you may need to manually provide your Matlab path to the script using -m option. Use -h option for more information)

---

- Otherwise, to get the outputs in Matlab, use the Matlab script with the command:
```objective
[lh_proj_data, rh_proj_data] = CBIG_RF_ProjectMNI2fsaverage('full/path/to/A.nii.gz');
```
  (Note that the outputs are not saved in this case)

---

Note that RF-ANTs MNI152-to-fsaverage mappings are used by default.
For more options (e.g. project data from Colin27 space, or use RF-FS mapping), please read the help message of the scripts.

  


