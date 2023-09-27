function [test_array] = CBIG_DiffProc_AMICO_check_example_results(test_struct)
% This function checks whether the example AMICO results are consistent
% with the generated example.
%
% Input:
% test_struct - A struct containing OD_dir and ICVF_dir, which refer to
%               the paths to the corresponding nifti files generated from AMICO.
%               E.g. s.OD_dir = fullfile(subj_output_dir, 'FIT_OD.nii.gz');
%                    s.ICVF_dir = fullfile(subj_output_dir, 'FIT_ICVF.nii.gz');
% Output:
% test_array - A array with 4 entries. 1 indicates a pass and 0 is a fail.
%              test_array(1): OD image volume sizes are equal
%              test_array(2): OD image values are similar
%              test_array(3): ICVF image volume sizes are equal
%              test_array(4): ICVF image values are similar
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% set directories
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
ref_dir = fullfile(CBIG_CODE_DIR,'stable_projects','preprocessing', ...
    'CBIG2022_DiffProc', 'examples', 'AMICO_example_data', 'ref_output', ...
    'output', 'AMICO', 'subject_01', 'AMICO', 'NODDI');

% compare results
test_array = zeros(1,4);
% read tbss skeletonised images and compare to reference
test_ODimage = MRIread(test_struct.OD_dir);
ref_ODimage = MRIread(fullfile(ref_dir, 'FIT_OD.nii.gz'));
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
ref_ICVFimage = MRIread(fullfile(ref_dir, 'FIT_ICVF.nii.gz'));
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