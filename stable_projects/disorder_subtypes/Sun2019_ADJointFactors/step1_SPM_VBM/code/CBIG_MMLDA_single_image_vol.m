function Vols = CBIG_MMLDA_single_image_vol(gm_image)
% Calculate the volume of single image.
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

V = spm_vol(gm_image);
Vols = zeros(numel(V),1);
for j=1:numel(V)
    tot = 0;
    for i=1:V(j).dim(3)
        img = spm_slice_vol(V(j),spm_matrix([0 0 i]),V(j).dim(1:2),0);
        img = img(isfinite(img));
        tot = tot + sum(img(:));
    end
    voxvol = abs(det(V(j).mat));
    Vols(j) = tot*voxvol;
end
end