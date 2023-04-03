Enclosed are the surface information based on Freesurfer's fsaverage6 surface mesh.

Some of the enclosed files are pre-generated as follows:
- `./surf/<lh/rh>.inflated.H`, `./surf/<lh/rh>.inflated.H`: you can generate these curvature files with `mris_curvature -w <lh/rh>.inflated`
- `./surf/<lh/rh>.sphere.left_right`: you can generate these curvature files with `mris_left_right_register -dist 1 lh.sphere rh.sphere lh.sphere.left_right rh.sphere.left_right`

Note that [`mris_left_right_register`](https://github.com/freesurfer/freesurfer/tree/dev/mris_left_right_register) is only available in Freesurfer 6.0.0 and above.