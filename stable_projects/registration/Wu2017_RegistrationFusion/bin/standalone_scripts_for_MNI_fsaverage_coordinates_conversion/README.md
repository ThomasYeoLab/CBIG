## Reference

Jianxiao Wu, Gia H. Ngo, Alexander Schaefer, Douglas Greve, Jingwei Li, Tong He, Bruce Fischl, Simon B. Eickhoff, B.T. Thomas Yeo. [**Accurate Nonlinear Mapping between MNI152/Colin27 Volumetric and FreeSurfer Surface Coordinate Systems**](http://people.csail.mit.edu/ythomas/publications/2018VolSurfMapping-HBM.pdf), *Human Brain Mapping*, 2018.

---

## Code Release

This folder contains standalone codes for using the Registration Fusion mappings to convert between MNI152 voxels and fsaverage vertices. The mappings are produced in the same way as those used in our paper, but using all 1490 GSP subjects.

To download the codes for replication of our paper, or for creating customised mappings, you would need to install our Github repository (Refer to https://github.com/ThomasYeoLab/CBIG). You can find the project under stable_projects/registration/Wu2017_RegistrationFusion.

---

## Before you start

Make sure that FreeSurfer and Matlab are already installed. In addition, please make sure that you have set up FreeSurfer properly (see https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall#Setup.26Configuration for instructions).

 To test if the scripts work, you can try the example (see Example section).

---

## Example

Run `CBIG_RF_conversion_example.sh` for an example. This converts an example set of RAS coordinates ([60, 0, 10]) in MNI152 space to the corresponding fsaverage vertex, and back to coordinates in MNI152 space.

(Note that you may need to manually provide your Matlab path to the script using -m option. Use -h option for more information)

The help messages of `CBIG_RF_MNICoord2fsaverageVertex.m` and `CBIG_RF_fsaverageVertex2MNICoord.m` also include examples for using each script alone.

--- 

## MNI152 coordinates to fsaverage vertices

To convert RAS coordinates ([x, y, z]) in MNI152 space to fsaverage vertex numbers, use the Matlab script with the command:

```objective
vertices = CBIG_RF_MNICoord2fsaverageVertex([x; y; z]);
```

To convert to coordinates in fsaverage surface space, the surface mesh should be specified (except for white surface mesh, which is used by default). For example, to get coordinates on the pial surface, use the following command:

```objective
[~, fs_coords]= CBIG_RF_MNICoord2fsaverageVertex([x; y; z], 'pial');
```

For multiple sets of coordinates, the input should be a matrix of size 3 x N, where N is the number of sets of coordinates. For example, change the input to `[x1, x2, x3; y1, y2, y3; z1, z2, z3]` for 3 sets of coordinates [x1, y1, z1], [x2, y2, z2] and [x3, y3, z3].

---

## fsaverage vertices to MNI152 coordinates

To convert fsaverage vertex numbers (e.g. vertex number x in left hemisphere) to RAS coordinates in MNI152 space, use the Matlab script with the command:

```objective
mni_coords = CBIG_RF_fsaverageVertex2MNICoord('lh', x);
```

For multiple vertex numbers, change the input x to an array of the vertex numbers. For example, use `1:163842` to get the corresponding coordinates for all the vertices in the left hemisphere.

For vertices in the right hemisphere, change the first argument to `'rh'`.

---

Note that the conversion of MNI152 coordiantes to fsaverage vertices uses fsaverage-to-MNI152 mapping, while the reverse conversion uses MNI152-to-fsaverage mapping. Therefore the results obtained may not be consistent with those from the final projection scripts (or the standalone scripts) to project data between MNI152 space and fsaverage surface.

  


