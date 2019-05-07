% Normalize PET with respect to cerebellum grey matter.
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

mri = MRIread(Mask_file);
vol = mri.vol;
mri_size = size(vol);
mask_vol1d = reshape(vol, [mri_size(1)*mri_size(2)*mri_size(3) 1]);

mri = MRIread(MRI_file);
vol = mri.vol;
mri_size = size(vol);
mri_vol1d = reshape(vol, [mri_size(1)*mri_size(2)*mri_size(3) 1]);

% average the cerebellum gm value in PET_in_T1
avg_cere_gm = mean(mri_vol1d(logical(mask_vol1d)));

% normalize the whole image by the average value
mri.vol = vol/avg_cere_gm;

MRIwrite(mri, Output_file);



