function [test_array] = CBIG_DiffProc_TBSS_check_example_result(test_struct)
% This function checks whether the example TBSS results are consistent for the
% unit test.
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects','preprocessing', ...
    'CBIG2022_DiffProc','TBSS');
% compare results
test_array = zeros(1,8);
% read tbss skeletonised images and compare to reference
test_FAimage = MRIread(test_struct.FA_dir);
ref_FAimage = MRIread(fullfile(ref_dir, 'ref_output', 'all_FA_skeletonised.nii.gz'));
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
ref_MDimage = MRIread(fullfile(ref_dir, 'ref_output', 'all_MD_skeletonised.nii.gz'));
% check size
if isequal(size(test_MDimage.vol), size(ref_MDimage.vol))
        test_array(5) = 1;
    else
        fprintf('Size of test and reference MD images are different.\n');
end
% check values for each subject
for n = 1:size(ref_MDimage.vol,4)
    if (abs(sum(test_MDimage.vol(:,:,:,n),'all') - sum(ref_MDimage.vol(:,:,:,n),'all')) < 1e-6)
            test_array(5+n) = 1;
        else
            fprintf('Test and reference MD values are different.\n');
    end
end

end