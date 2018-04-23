function CBIG_RF_make_cortexMask(template_sub_id, count_map, num_sub, output_dir)
% CBIG_RF_make_cortexMask(template_sub_id, count_map, num_sub, output_dir)
%
% This function creates loose cortical mask for a volumetric atlas, 
% using surface-to-volume count maps
% The mask is created by thresholding the mapping by 15% of total subject used
%
% Input:
%     - template_sub_id:
%                        recon-all subject ID of the volumetric atlas
%     - count_map      :
%                        surface-to-volume count map created when computing 
%                          surface-to-volume average map
%     - num_sub        :
%                        number of subjects used when computing the average map
%     - output_dir     :
%                        absolute path to directory where output should be stored
%
% Output:
%     - There is no function output.
%     - 1 volume mask is created in output_dir:
%           [output_prefix].nii.gz
%
% Example:
% CBIG_RF_make_cortexMask('FSL_MNI152_FS4.5.0, 'results/allSub_fsaverage_to_FSL_MNI152_FS4.5.9_RF_ANTs_count.mat, 
%                       1490, 'results')
% This command generates cortical mask for the FSL_MNI152_FS4.5.0 template using RF-ANTs projection 
% count map across 1490 subject
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin ~= 4
    disp('usage: CBIG_RF_make_cortexMask(template_sub_id, count_map, num_sub, output_dir)');
    return
end

%Read Freesurfer aparc+aseg parcellation
aparc_aseg = MRIread([getenv('CBIG_CODE_DIR') '/data/templates/volume/' ...
    template_sub_id '/mri/aparc+aseg.mgz']);
%Compute cortical and noncortical masks from aparc+aseg
mask_fs = (aparc_aseg.vol >= 1000);
mask_nc = ((aparc_aseg.vol < 1000) & (aparc_aseg.vol > 0) & ...
    (aparc_aseg.vol ~= 41) & (aparc_aseg.vol ~= 2));

%Read and threshold count map
thr = num_sub * 0.15; %theshold at 15%
load(count_map);
mask_count = (lh_count >= thr) + (rh_count >= thr);
mask_count = reshape(mask_count, size(aparc_aseg.vol));

%Merge masks and smooth
mask = logical(mask_fs) | logical(mask_count);
mask = smooth3(double(mask), 'box', 5);
mask(mask < 0.5) = 0;
mask(mask >= 0.5) = 1;
mask = logical(mask);

%Remove non-cortical regions
mask(mask_nc) = 0;

%Fill holes and remove islands
mask = imfill(mask, 'holes');
[l, num] = bwlabeln(mask);
for i = 1:num
    if (sum(l(:)==i) < 100)
        mask(l==i) = 0;
    end
end

%Save results
output = aparc_aseg;
output.vol = mask;
MRIwrite(output, [output_dir '/' template_sub_id '_cortex_estimate.nii.gz']);

end
