function [test_array] = CBIG_DiffProc_AMICO_check_example_result(test_struct)
% This function checks whether the example AMICO results are consistent for the
% unit test.
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects','preprocessing', ...
    'CBIG2022_DiffProc','AMICO');
% compare results
test_array = zeros(1,4);
% read tbss skeletonised images and compare to reference
test_ODimage = MRIread(test_struct.OD_dir);
ref_ODimage = MRIread(fullfile(ref_dir, 'ref_output', 'FIT_OD.nii.gz'));
% check size
if isequal(size(test_ODimage.vol), size(ref_ODimage.vol))
        test_array(1) = 1;
    else
        fprintf('Size of test and reference OD images are different.\n');
end
% check values for each subject
if (abs(sum(test_ODimage.vol,'all') - sum(ref_ODimage.vol,'all')) < 1e-6)
        test_array(2) = 1;
    else
        fprintf('Test and reference OD values are different.\n');
end

% read tbss skeletonised images and compare to reference
test_ICVFimage = MRIread(test_struct.ICVF_dir);
ref_ICVFimage = MRIread(fullfile(ref_dir, 'ref_output', 'FIT_ICVF.nii.gz'));
% check size
if isequal(size(test_ICVFimage.vol), size(ref_ICVFimage.vol))
        test_array(3) = 1;
    else
        fprintf('Size of test and reference ICVF images are different.\n');
end
% check values for each subject
if (abs(sum(test_ICVFimage.vol,'all') - sum(ref_ICVFimage.vol,'all')) < 1e-6)
        test_array(4) = 1;
    else
        fprintf('Test and reference ICVF values are different.\n');
end

end