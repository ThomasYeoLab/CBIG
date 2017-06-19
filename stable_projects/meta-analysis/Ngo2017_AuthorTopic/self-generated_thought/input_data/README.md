`MNI152_activation_coordinates` contains the activation coordinates in MNI152 space of the 7 task involving self-generated thought.

`self-generated_thought.mat` is the Matlab input used by the inference algorithm. Each experimental contrast corresponds to a 2-mm-resolution binary activation image, in which a voxel was given a value of 1 if it was within 10mm of any activation focus, and 0 otherwise. See Appendix A5 in the paper for more details.

The functions used to produce the final Matlab input from the original activation coordinates will be released soon.
