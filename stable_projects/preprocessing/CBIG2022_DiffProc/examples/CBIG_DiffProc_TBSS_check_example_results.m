function [test_array] = CBIG_DiffProc_TBSS_check_example_results(test_struct)
% This function checks whether the example TBSS results are consistent
% with the generated example.
%
% Input:
% test_struct - A struct containing FA_dir and MD_dir, which refer
%               to the paths to the corresponding FA and MD TBSS skeletons.
%               E.g. s.FA_dir = fullfile('stats', 'all_FA_skeletonised.nii.gz');
%                    s.MD_dir = fullfile('stats', 'all_MD_skeletonised.nii.gz');
%
% Output:
% test_array - A array with 5 entries. 1 indicates a pass and 0 is a fail.
%              test_array(1): size of FA skeleton are equal
%              test_array(2): FA values for subject 1 are similar
%              test_array(3): FA values for subject 2 are similar
%              test_array(4): size of MD skeleton are equal
%              test_array(5): MD values for subject 1 are similar
%              test_array(6): MD values for subject 2 are similar
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% set directories
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
ref_dir = fullfile(CBIG_CODE_DIR,'stable_projects', 'preprocessing', ...
        'CBIG2022_DiffProc', 'examples', 'TBSS_example_data', ...
        'ref_output');
    
% compare results
test_array = zeros(1,6);
% read tbss skeletonised images and compare to reference
test_FAimage = MRIread(test_struct.FA_dir);
ref_FAimage = MRIread(fullfile(ref_dir, 'stats', 'all_FA_skeletonised.nii.gz'));
% check size
if isequal(size(test_FAimage.vol), size(ref_FAimage.vol))
        test_array(1) = 1;
    else
        fprintf('Size of test and reference FA images are different.\n');
end
% check values for each subject
for n = 1:size(ref_FAimage.vol,4)
    if (abs(sum(test_FAimage.vol(:,:,:,n),'all') - sum(ref_FAimage.vol(:,:,:,n),'all')) < 1e-6)
            test_array(1+n) = 1;
        else
            fprintf('Test and reference FA values are different.\n');
    end
end

% clear images for space
clear test_FAimage ref_FAimage

% read tbss skeletonised images and compare to reference
test_MDimage = MRIread(test_struct.MD_dir);
ref_MDimage = MRIread(fullfile(ref_dir, 'stats', 'all_MD_skeletonised.nii.gz'));
% check size
if isequal(size(test_MDimage.vol), size(ref_MDimage.vol))
        test_array(4) = 1;
    else
        fprintf('Size of test and reference MD images are different.\n');
end
% check values for each subject
for n = 1:size(ref_MDimage.vol,4)
    if (abs(sum(test_MDimage.vol(:,:,:,n),'all') - sum(ref_MDimage.vol(:,:,:,n),'all')) < 1e-6)
            test_array(4+n) = 1;
        else
            fprintf('Test and reference MD values are different.\n');
    end
end

end