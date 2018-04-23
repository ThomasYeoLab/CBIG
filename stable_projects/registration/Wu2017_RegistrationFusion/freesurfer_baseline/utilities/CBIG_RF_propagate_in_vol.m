function [output, output_seg] = CBIG_RF_propagate_in_vol(input_lh, input_rh, mask_input)
% CBIG_RF_propagate_in_vol(input_lh, input_rh, mask_input)
%
% This function grows a volume (e.g. segmentation labels) to fill up an input mask
%
% Input:
%     - input_lh  :
%                   path to input volume in left hemisphere
%     - input_rh  :
%                   path to input volume in right hemisphere
%     - mask_input:
%                   path to the input mask which defines regions to be filled 
%                   (e.g. a liberal cortical mask)
%                   (default: $CBIG_CODE_DIR/bin/liberal_cortex_masks/
%                     MNI152_norm_cortex_estimate.nii.gz)
%
% Output:
%     - output    :
%                   output volume after resolving midline issues 
%                   (no propagation)
%     - output_seg:
%                   propagated output volume in segmentation format and masked
%                   (rh values have a base of 1000)
%
% Example:
% [output, output_seg] = CBIG_RF_propagate_in_vol('lh.fsaverage2MNI_annot.nii.gz', 
%                 'rh.fsaverage2MNI_annot.nii.gz', 'MNI_cortex_estimate.nii.gz')
% This command takes in the two input volumes and returns the output by removing data in the midline voxels. 
% The volume will then be propagated to fill the MNI cortex mask. 
% The sementation formatted volume output of this propagated volume is then returned as output_seg.
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%Function usage
if nargin < 2
    disp('usage: CBIG_RF_propagate_in_vol(input_lh, input_rh, mask_input');
    return
end

%Default parameter
dir_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 3
    mask_input = [dir_root '/bin/liberal_cortex_masks/MNI152_norm_cortex_estimate.nii.gz'];
end

%Read inputs
lh = MRIread(input_lh);
rh = MRIread(input_rh);
mask = MRIread(mask_input);

%Resolve midline issues
midline = double(lh.vol==0) + double(rh.vol==0);
lh.vol(midline==0) = 0;
rh.vol(midline==0) = 0;

%Construct output
output = lh;
output.vol = lh.vol + rh.vol;

%Construct segmentation formatted output
output_seg = output;
output_seg.vol(rh.vol~=0) = output_seg.vol(rh.vol~=0) + 1000;
[~, id] = bwdist(output_seg.vol);
output_seg.vol = output_seg.vol(id);
output_seg.vol(mask.vol==0) = 0;

end
