function fs5 = CBIG_pFIC_ROI2fs5(roi, roi_list)

% fs5 = CBIG_pFIC_ROI2fs5(roi, roi_list)
% This funciton upsamples the ROI-level values to fsaverage5 vertex-level values.
% roi_list is to specify which ROI(s) are to exclude from upsampling.
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Note that this function only upsamples from Deiskan parcellation to fsaverage5 space
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% Input:
%   - roi: a 68-by-1 vector according to Desikan parcellation
%   - roi_list: a 68-by-1 binary vector. 1 means the ROI is to be upsampled to
%   fsaverage5 as per normal; 0 means the ROI is to be replaced by 0 during
%   upsampling, this is typically used when certain ROIs are excluded due
%   to failing a statistical test
% Output:
%   - fs5: a 20484-by-1 vector (left hemisphere followed by right
%   hemisphere). This is an upsampled vertex-level version of the input
%   roi-level input values.
% Example:
% EI_contrast_fs5 = CBIG_pFIC_ROI2fs5(EI_contrast_roi, [ones(31,1); 0; ones(36, 1)]);
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

disp('Upsampling from Desikan parcellation to fsaverage5 space...')
lh_cortex = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
lh_cortex = lh_cortex.MARS_label;
rh_cortex = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'cortex');
rh_cortex = rh_cortex.MARS_label;
whole_cortex = [lh_cortex rh_cortex]';

fs5 = nan(20484, 1);

if length(roi_list) == 68
    roi(~roi_list) = inf;
else
    disp('ROI list length is not right!')
    return
end
lh_desikan = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'aparc.annot');
lh_desikan = lh_desikan.MARS_label;
rh_desikan = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'aparc.annot');
rh_desikan = rh_desikan.MARS_label + 36;

whole_desikan = [lh_desikan rh_desikan]';

roi = [0; roi(1:3); 0; roi(4:34); 0; roi(35:37); 0; roi(38:68)];

for i = 1:20484
   fs5(i) = roi(whole_desikan(i)); 
   if whole_cortex(i) == 1
      fs5(i) = NaN; 
   end
end

fs5(isnan(fs5)) = [];
disp('Upsample to fsaverage5 done!')

end
