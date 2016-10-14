This folder contains small, but useful data files required by various functions in the repository, for example, MNI152 templates that have been preprocessed through Freesurfer's `recon-all` command.

Note that there are two symbolic links `surf2surf_gui_data`, and `vol2surf_gui_data` that point to data required by interactive GUI tools hosted on the server at Thomas' lab (not available to the public).
`surf2surf_gui_data` is required by a GUI tool (`Surf2SurfGui.m`) that interactively shows functional connectivity within the surface.
`vol2surf_gui_data` is required by a GUI tool (`Vol2SurfGui.m)` that interactively shows functional connectivity between volume and surface.
