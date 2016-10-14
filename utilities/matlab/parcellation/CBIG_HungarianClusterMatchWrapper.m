function CBIG_HungarianClusterMatchWrapper(ref_file, input_file, output_file)

% CBIG_HungarianClusterMatchWrapper(ref_file, input_file, output_file)
%
% This is the wrapper function of Hungarian matching for labeling.
%
% Given a reference file name, this function match labels of input file
% with the reference file.
%
% Reference file and input file are data in mgh, mgz, img, bhdr, nii, or
% nii.gz format.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



ref_labels = MRIread(ref_file);
input_labels = MRIread(input_file);

output = input_labels;
[tmp, assign, cost, dice_overlap] = CBIG_HungarianClusterMatch(ref_labels.vol(:), input_labels.vol(:));
output.vol = reshape(tmp, size(input_labels.vol));
disp(num2str(cost));
disp(num2str(dice_overlap));
MRIwrite(output, output_file);

exit;
